# Plot LP results for one week
# Meant to complement (and is almost entirely based off of) plot_sam_results.jl
# Also takes code from simpleouterloop.jl
# 2018-03-18
#include("module_innerloop.jl")
include("module_NLP_loop.jl")

# Simple outer loop
#using inner_loop
using NLP_loop
using DataFrames
import CSV
using Gadfly
using JuMP
using AmplNLWriter
using Clp


############################ OUTER LOOP ####################################

equip_data = readcsv("Equipment_Costs.csv")[2:end,:]
damper_data = readcsv("Damper_Data.csv")[2:end,:]

equip_costs = equip_data[:, 2]
equip_capacity = equip_data[:, 3]
equip_efficiency = equip_data[:, 4]
ISMRE = equip_data[:, 5]

damper_placement = damper_data[:, 2]
damper_cost = 5000

optimal_cost = Inf
best = 0



nYears = 20
discountRate = 0.05

DO1 = [1,1,1,1,1,1,1,1,1,1,1,1]
DO2 = [0,1,1,1,1,1,1,1,1,1,1,1]
DO3 = [0,0,1,1,1,1,1,1,1,1,1,1]
DO4 = [0,0,0,1,1,1,1,1,1,1,1,1]
DO5 = [0,0,0,0,1,1,1,1,1,1,1,1]
DO6 = [0,0,0,0,0,1,1,1,1,1,1,1]
DO7 = [0,0,0,0,0,0,1,1,1,1,1,1]
DO8 = [0,0,0,0,0,0,0,1,1,1,1,1]
DO9 = [0,0,0,0,0,0,0,0,1,1,1,1]
DO10 = [0,0,0,0,0,0,0,0,0,1,1,1]
DO11 = [0,0,0,0,0,0,0,0,0,0,1,1]
DO12 = [0,0,0,0,0,0,0,0,0,0,0,1]
DO13 = [0,0,0,0,0,0,0,0,0,0,0,0]

damperPlacementOptions = [DO1, DO2, DO3, DO4, DO5, DO6, DO7, DO8, DO9, DO10, DO11, DO12, DO13]

# Note that right now our inner loop is implemented such that the efficiency
# number it takes is in kWh/m^3 where more efficiency numbers are given in
# m^3/kWh, so we take the inverse of the efficiency when we pass it in on line
# 33.




#=
feasibility, opex, fmax, pmax, hmax, NLPstatus =
    OpsOptLP(T, equip_efficiency[i], equip_capacity[i], ISMRE[i], damperPlacementOptions[j])

capex = equip_costs[i] + sum(damper_cost * (1 - damperPlacementOptions[j][k]) for k = 1:length(damperPlacementOptions[j]))

NPC = capex + sum(opex/(1+discountRate)^p for p=1:nYears)
=#
########################### Inner Loop ########################################
PATH_TO_SOLVERS = ENV["ERE291_SOLVERS"]

#function OpsOptLP(T, P, capacity, ISMRE, damperPlacement)
i = 1
j = 13

equip_efficiency[i]
ISMRE[i]

T = 168
P = 0.00045 # kWh/m^3
capacity = equip_capacity[i]
ISMRE = 0.5
damperPlacement = damperPlacementOptions[j]



###############################
######### Define Sets #########
###############################

T_repeat = 24

data = readcsv("Room-Data_v3.0.csv")
AQdata = readcsv("AirQualityData2016.csv")
rooms = data[2:end,3]
N = length(rooms)

###################################################
############ Define parameters and data ###########
###################################################

### Air Quality Parameters ###

# CO2 concentration of outdoor air at time t [ppm]
AMB_CO2 = AQdata[3:end,3]

# PM2.5 concentration of outdoor air at time t [ug / m3]
AMB_PM2_5 = AQdata[3:end,2]

# moisture ratio of outdoor air at time t [gm water / gm of Dry Air]
AMB_HUMID = AQdata[3:end,4]

# maximum allowable indoor CO2 concentration [ppm]
CO2_MAX = 1000

# maximum allowable indoor PM2.5 concentration [ug/m3]
PM25_MAX = 35

# maximum allowable indoor humidity ratio [gm water / gm of Dry Air]
HUMID_MAX = 0.012

# Calculate rate of CO2 emissions from occupancy at time t [cm3 / hr]
massPerson = 150 #lbs
met = 0.843318 #kcal/hr-lb
qCO2 = 241 #cm3/kcal

CO2pp = massPerson * met * qCO2 #cm3 CO2 / person-hr

humidityPerPerson = 0.15 #[kg / person-hr]

roomData = Dict{Int64,Array{Float64}}()
roomOccupancy = Dict{Int64, Array{Int64}}() #number of People in each hour
roomCO2Source = Dict{Int64, Array{Float64}}() #cm3 of CO2 Produced in each hour
roomCO2Max = Dict{Int64, Float64}() #cm3 of CO2 Produced in each hour
roomHumidSource = Dict{Int64, Array{Float64}}() #kg H20 produced per hour

for i = 1:N
    roomData[rooms[i]] = data[i+1,4:7]
    roomOccupancy[rooms[i]] = data[i+1,8:31]
    roomCO2Source[rooms[i]] = CO2pp .* roomOccupancy[rooms[i]]
    roomHumidSource[rooms[i]] = humidityPerPerson .* roomOccupancy[rooms[i]]
end

CMHOrig = zeros(N,24)
CMHpp = 55 #CMH per person, US guideline of 15 cfm/person

for i = 1:N
    CMHOrig[i,:] = CMHpp .* roomOccupancy[rooms[i]]
end

CMHOptResult, fmax, NLPstatus = CMH_opt(T_repeat, N, P, damperPlacement)

CMH = zeros(N,T_repeat)

for i in 1:N
    for j in 1:T_repeat
        CMH[i,j] = max(CMHOrig[i,j], CMHOptResult[i,j])
    end
end

CMH = repmat(CMH, 1, ceil(Int, T / T_repeat))

# initial indoor CO2 concentration at t = 0
CO2_0 = AMB_CO2[1]

# initial indoor PM2.5 concentration at t = 0
PM25_0 = AMB_PM2_5[1]

# initial indoor PM2.5 concentration at t = 0
HUMID_0 = AMB_HUMID[1]

# density of air [kg / m3]
airDensity = 1.225

### Equipment Parameters ###

# PM2.5 absorption capacity per filter before replacement [ug/filter]
k = 10^7

# pressure drop through ducts [Pa]
pressureDropDucts = 103

# pressure drop through 1x HEPA filter [{Pa]
pressureDropHEPAFilter = 250

# pressure drop total [Pa]
pressureDropSystem = pressureDropDucts + pressureDropHEPAFilter


### Cost Parameters ###

# cost of electricity at time t [RMB/kWh]
Celec = 1.4

# cost per in-room PM2.5 filter [RMB/filter]
Cfilter = 50


####################################
######### Initialize Model #########
####################################

#m = Model(solver=AmplNLSolver(joinpath(PATH_TO_SOLVERS,"knitro"), ["outlev=0", "ms_enable=1"]))
m = Model(solver = ClpSolver())

####################################
######## Decision variables ########
####################################

# rate of PM2.5 absorption in room i at time t [ug / hr]
@variable(m, roomPM25Absorption[1:N, 1:T] >= 0)

# rate of humidity removal at time t [kg / hr]
@variable(m, kgMoistureRemoved[1:N, 1:T] >= 0)

# CO2 concentration in room i at time t [m3 / ug]
@variable(m, CO2[1:N, 1:T] >= 0)

# PM2.5 concentration in room i at time t [ppm]
@variable(m, PM25[1:N, 1:T] >= 0)

# humidity ratio of indoor air at time t [gm water / gm of dry air]
@variable(m, HUMID[1:N, 1:T] >= 0)


#####################################
######## Dependent variables ########
#####################################

# number of PM2.5 filters required for all rooms
@expression(m, numFilters, sum(roomPM25Absorption) / k)

# rate of incoming humidity from FAU at time t [kg / hr]
@expression(m, kgMoistureFAUIn[i=1:N, t=1:T], CMH[i, (t-1)%24 + 1] * airDensity * AMB_HUMID[t])

# rate of outgoing humidity from FAU at time t [kg / hr]
@expression(m, kgMoistureFAUOut[i=1:N, t=1:T], CMH[i, (t-1)%24 + 1] * airDensity * HUMID[i,t])

######################################
######## Objective Functions #########
######################################

# Minimize the total cost, equal to the sum of the cost of fan electricity over the operational period, plus the cost of PM2.5 filter replacements
@objective(m, Min, Celec*((P + pressureDropSystem / 3600 / 1000)*sum(CMH) + ISMRE*sum(kgMoistureRemoved)) + Cfilter*numFilters)


######################################
############# Constraints ############
######################################

# CO2 concentration in the room at time t is equal to the CO2 concentration in the room at t-1 plus (the CO2 mass introduced at t minus the CO2 mass removed at t) divided by the room volume
# The constraint is structured this way because the HVAC system can only react to CO2 levels at time t, it cannot pre-emptively condition the air at t-1
# What if we just assumed an air change rate based on the maximum CO2
@constraint(m, [i=1:N, t=2:T], CO2[i,t] == CO2[i,t-1] + (roomCO2Source[rooms[i]][(t-1)%24 + 1] + CMH[i,(t-1)%24 + 1]*(AMB_CO2[t] - CO2[i,t]))/roomData[rooms[i]][2])

# PM2.5 concentration in the room at time t is equal to the PM 2.5 in the room at t-1 plus the mass of PM2.5 introduced at time t minus the mass filtered at time t divided by the room volume
@constraint(m, [i=1:N, t=2:T], PM25[i,t] == PM25[i,t-1] + (CMH[i, (t-1)%24 + 1]*(AMB_PM2_5[t] - PM25[i,t]) - roomPM25Absorption[i,t])/roomData[rooms[i]][2])

# Humidity ratio in the room at time t is equal to the Humidity ratio  in the room at t-1 plus (the H2O mass introduced at t minus the H2O mass removed at t) divided by the room volume
# The constraint is structured this way because the HVAC system can only react to humidity levels at time t, it cannot pre-emptively condition the air at t-1
@constraint(m, [i=1:N, t=2:T], HUMID[i,t] == HUMID[i,t-1] + (roomHumidSource[rooms[i]][(t-1)%24 + 1] + kgMoistureFAUIn[i,t] - kgMoistureFAUOut[i,t] - kgMoistureRemoved[i,t])/(roomData[rooms[i]][2] * airDensity))

# set  initial condition to ambient concentrations
@constraint(m, [i=1:N], PM25[i,1] == PM25_0)
@constraint(m, [i=1:N], CO2[i,1] == CO2_0)
@constraint(m, [i=1:N], HUMID[i,1] == HUMID_0)

# CO2, PM2.5 and Humidity concentrations reset at the end of the cycle.
# @constraint(m, [i=1:N], PM25[i,T] == PM25[i,1])
# @constraint(m, [i=1:N], CO2[i,T] == CO2[i,1])
# @constraint(m, [i=1:N], HUMID[i,T] == HUMID[i,1])

# maximum allowable indoor CO2 concentration [ppm].  Constraints do not apply to first timestep, since conditioning has not yet been applied.
@constraint(m, [i=1:N, t=2:T], CO2[i,t] <= CO2_MAX)

# maximum allowable indoor PM2.5 concentration [ug/m3].  Constraints do not apply to first timestep, since conditioning has not yet been applied.
@constraint(m, [i=1:N, t=2:T], PM25[i,t] <= PM25_MAX)

# maximum allowable indoor humidity ratio [gm water / gm of Dry Air].  Constraints do not apply to first timestep, since conditioning has not yet been applied.
@constraint(m, [i=1:N, t=2:T], HUMID[i,t] <= HUMID_MAX)

######################################
########### Print and solve ##########
######################################

#print(m)
solve(m)





CO2output = getvalue(CO2)
PM25output = getvalue(PM25)
HUMIDoutput = getvalue(HUMID)

CO2result = zeros(1,T)
PM25result = zeros(1,T)
HUMIDresult = zeros(1,T)

for t = 1:T
    CO2result[t] = sum(CO2output[i,t] * roomData[rooms[i]][2] for i in 1:N) / sum(roomData[rooms[i]][2] for i in 1:N)
    PM25result[t] = sum(PM25output[i,t] * roomData[rooms[i]][2] for i in 1:N) / sum(roomData[rooms[i]][2] for i in 1:N)
    HUMIDresult[t] = sum(HUMIDoutput[i,t] * roomData[rooms[i]][2] for i in 1:N) / sum(roomData[rooms[i]][2] for i in 1:N)
end


PM25Absorbedresult = sum(getvalue(roomPM25Absorption),1)
HUMIDAbsorbedresult = sum(getvalue(kgMoistureRemoved),1)
CMH_total = sum(CMH, 1)

LP_data = DataFrame(
    Timestep = 1:T,
    CMH = convert(DataFrame, CMH_total')[:x1],
    CO2 = convert(DataFrame, CO2result')[:x1],
    PM25_removed = convert(DataFrame, PM25Absorbedresult')[:x1],
    PM25 = convert(DataFrame, PM25result')[:x1],
    Humidity_removed = convert(DataFrame, HUMIDAbsorbedresult')[:x1],
    Humidity = convert(DataFrame, HUMIDresult')[:x1]
    )


########################## Plotting ###########################################

light_theme = Theme(
    background_color = "white",
    panel_fill = "white",
    default_color = "blue",
    key_position = :right
    )

CMH_plot =
    plot(LP_data,
        x = "Timestep",
        y = "CMH",
        Geom.line,
        Guide.Title("Air Flow into Building"),
        Guide.XLabel("Time (hrs)"),
        Guide.YLabel("Air Flow (CMH)"),
        light_theme
    )

CO2_plot =
    plot(
        LP_data,
        x = "Timestep",
        y = "CO2",
        Geom.line,
        Guide.Title("CO2 Concentration inside room (ppm)"),
        Guide.XLabel("Time (hrs)"),
        Guide.YLabel("CO2 Concentration (ppm)"),
        light_theme
    )

PM25_plot =
    plot(
        LP_data,
        x = "Timestep",
        y = "PM25",
        Geom.line,
        Guide.Title("PM2.5 Concentration (µg)"),
        Guide.XLabel("Time (hrs)"),
        Guide.YLabel("PM2.5 Concentration (µg / m^3)"),
        light_theme
    )

PM25_absorbed_plot =
    plot(
        LP_data,
        x = "Timestep",
        y = "PM25_removed",
        Geom.line,
        Guide.Title("PM2.5 Removed by Filters"),
        Guide.XLabel("Time (hrs)"),
        Guide.YLabel("PM2.5 Removed (µg)"),
        light_theme
    )

HUMID_plot =
    plot(
        LP_data,
        x = "Timestep",
        y = "Humidity",
        Geom.line,
        Guide.Title("Humidity"),
        Guide.XLabel("Time (hrs)"),
        Guide.YLabel("Humidity Ratio (moisture content)"),
        light_theme
    )

HUMID_absorbed_plot =
    plot(
        LP_data,
        x = "Timestep",
        y = "Humidity_removed",
        Geom.line,
        Guide.Title("Moisture Removed by Dehumidification"),
        Guide.XLabel("Time (hrs)"),
        Guide.YLabel("Moisture Removed (g H2O)"),
        light_theme
    )


final = vstack(
    hstack(CMH_plot, CO2_plot),
    hstack(PM25_plot, HUMID_plot)
    #hstack(HUMID_plot, HUMID_absorbed_plot)
    )

img = PNG("LP_one_week_short.png", 12inch, 12inch, dpi = 1200)
draw(img, final)
