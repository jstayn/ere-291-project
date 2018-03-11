####################################
######### Initialize tools #########
####################################

PATH_TO_SOLVERS = ENV["ERE291_SOLVERS"]

using JuMP
using Clp
using Cbc
using AmplNLWriter

using DataFrames
using Gadfly


###############################
######### Define Sets #########
###############################

# Currently testing a period of one year (8760) hours
T = 8760


###################################################
############ Define parameters and data ###########
###################################################

### Air Quality Parameters ###

# Calculate rate of CO2 emissions from occupancy at time t [cm3 / hr]

data = readcsv("Master-Data_v2.0.csv")
AQdata = readcsv("AirQualityData2016.csv")
rooms = data[2:end,3]

# CO2 concentration of outdoor air at time t [ppm]
AMB_CO2 = AQdata[3:end,3]

# PM2.5 concentration of outdoor air at time t [ug / m3]
AMB_PM2_5 = AQdata[3:end,2]

# moisture ratio of outdoor air at time t [gm water / gm of Dry Air]
AMB_HUMID = AQdata[3:end,4]

# maximum allowable indoor CO2 concentration [ppm]
CO2_MAX = 800

# maximum allowable indoor PM2.5 concentration [ug/m3]
PM25_MAX = 35

# maximum allowable indoor humidity ratio [gm water / gm of Dry Air]
HUMID_MAX = 0.012

massPerson = 150 #lbs
met = 0.843318 #kcal/hr-lb
qCO2 = 241 #cm3/kcal

CO2pp = massPerson * met * qCO2 #cm3 CO2 / person-hr

roomData = Dict{Int64,Array{Float64}}()
roomOccupancy = Dict{Int64, Array{Int64}}() #number of People in each hour
roomCO2Source = Dict{Int64, Array{Float64}}() #cm3 of CO2 Produced in each hour
roomCO2Max = Dict{Int64, Float64}() #cm3 of CO2 Produced in each hour
CMH = zeros(24)

for i = 1:length(rooms)
    roomData[rooms[i]] = data[i+1,4:6]
    roomOccupancy[rooms[i]] = data[i+1,7:30]
    roomCO2Source[rooms[i]] = CO2pp .* roomOccupancy[rooms[i]]
    roomCO2Max[rooms[i]] = maximum(roomCO2Source[rooms[i]])
end

# Set CMH during occupied hours based on the difference between the ppm generated and the max allowable ppm level

### INSERT FUNCTION CALL TO 24-HOUR NONLINEAR CMH SOLVER HERE ###
### SHOULD RETURN CMH PER ROOM FOR A SINGLE DAY ###
### REPLACE HARD-CODED CMH PARAMETER ###

ACH = 5.15
CMH = [0, 0, 0, 0, 0, 0, 0, 0, 0, roomData[101][2]*ACH, roomData[101][2]*ACH, roomData[101][2]*ACH, roomData[101][2]*ACH, roomData[101][2]*ACH,
    roomData[101][2]*ACH, roomData[101][2]*ACH, roomData[101][2]*ACH, roomData[101][2]*ACH, roomData[101][2]*ACH, 0, 0, 0, 0, 0, 0]


# Calculate rate of humidity generated from occupancy at time t [kg / hr]

humidityPerPerson = 0.20 #[kg / hr]
roomHumidSource = Dict{Int64, Array{Float64}}()

for i = 1:length(rooms)
    roomData[rooms[i]] = data[i+1,4:6]
    roomOccupancy[rooms[i]] = data[i+1,7:30]
    roomHumidSource[rooms[i]] = humidityPerPerson .* roomOccupancy[rooms[i]]
end

println(roomHumidSource[101])

# initial indoor CO2 concentration at t = 0
CO2_0 = AMB_CO2[1] #do we still need this? - YP

# initial indoor PM2.5 concentration at t = 0
PM25_0 = AMB_PM2_5[1] #do we still need this? - YP

# initial indoor PM2.5 concentration at t = 0
HUMID_0 = AMB_HUMID[1] #do we still need this? - YP

# density of air [kg / m3]
airDensity = 1.225

### Equipment Parameters ###

# max intake capacity of FAU model i [m3 air / hr]
x = 2000 #dummy varaible - need to merge with John's design loop - YP

# efficiency of the fan of FAU model i [kWh / m3 air]
P = .00045

# PM2.5 absorption capacity per filter before replacement [ug/filter]
k = 10^7

# minimum integrated seasonal moisture removal efficiency (ISMRE) [kWh/kg]
ISMRE = 0.5

# efficiency of PM2.5 filter
R = 0.8 #do we still need this? - YP


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

# rate of air intake at time t [m3 / hr]
# @variable(m, CMH[1:T] >= 0, start=1) #need to change to absorption per room - YP

# rate of PM2.5 absorption in room i at time t [ug / hr]
@variable(m, roomPM25Absorption[1:T] >= 0)

# rate of humidity removal at time t [kg / hr]
@variable(m, kgMoistureRemoved[1:T] >= 0)

# CO2 concentration in room i at time t [m3 / ug]
@variable(m, CO2[1:T] >= 0) #need to change to CO2concentration per room - YP

# PM2.5 concentration in room i at time t [ppm]
@variable(m, PM25[1:T] >= 0) #need to change to PM2.5 concentration per room - YP

# humidity ratio of indoor air at time t [gm water / gm of dry air]
@variable(m, HUMID[1:T] >= 0) #need to change to PM2.5 concentration per room - YP


#####################################
######## Dependent variables ########
#####################################

# number of PM2.5 filters required for room i at time t [kg / hr]
@expression(m, N, 1/k*sum(roomPM25Absorption[t] for t in 1:T)) #need to change to filters per room - YP

# rate of incoming humidity from FAU at time t [kg / hr]
@expression(m, kgMoistureFAUIn[t=1:T], CMH[t%24 + 1] * airDensity * AMB_HUMID[t])

# rate of outgoing humidity from FAU at time t [kg / hr]
@expression(m, kgMoistureFAUOut[t=1:T], CMH[t%24 + 1] * airDensity * HUMID[t])

######################################
######## Objective Functions #########
######################################

# Minimize the total cost, equal to the sum of the cost of fan electricity over the operational period, plus the cost of PM2.5 filter replacements
@objective(m, Min, Celec*sum(CMH[t%24 + 1]*P + kgMoistureRemoved[t]*ISMRE for t in 1:T) + Cfilter*N)


######################################
############# Constraints ############
######################################

# CO2 concentration in the room at time t is equal to the CO2 concentration in the room at t-1 plus (the CO2 mass introduced at t minus the CO2 mass removed at t) divided by the room volume
# The constraint is structured this way because the HVAC system can only react to CO2 levels at time t, it cannot pre-emptively condition the air at t-1
# What if we just assumed an air change rate based on the maximum CO2
@constraint(m, [t=2:T], CO2[t] == CO2[t-1] + (roomCO2Source[101][t%24 + 1] + CMH[t%24 + 1]*AMB_CO2[t] - CMH[t%24 + 1]*CO2[t])/roomData[101][2])

# PM2.5 concentration in the room at time t is equal to the PM 2.5 in the room at t-1 plus the mass of PM2.5 introduced at time t minus the mass filtered at time t divided by the room volume
@constraint(m, [t=2:T], PM25[t] == PM25[t-1] + (CMH[t%24 + 1]*AMB_PM2_5[t] - CMH[t%24 + 1]*PM25[t] - roomPM25Absorption[t])/roomData[101][2])

# Humidity ratio in the room at time t is equal to the Humidity ratio  in the room at t-1 plus (the H2O mass introduced at t minus the H2O mass removed at t) divided by the room volume
# The constraint is structured this way because the HVAC system can only react to humidity levels at time t, it cannot pre-emptively condition the air at t-1
@constraint(m, [t=2:T], HUMID[t] == HUMID[t-1] + (roomHumidSource[101][t%24 + 1] + kgMoistureFAUIn[t] - kgMoistureFAUOut[t] - kgMoistureRemoved[t])/(roomData[101][2] * airDensity))

# set dummy initial condition to ambient concentrations
# @constraint(m, PM25[1] == PM25_0)
# @constraint(m, CO2[1] == CO2_0)
# @constraint(m, HUMID[1] == HUMID_0)

# CO2, PM2.5 and Humidity concentrations reset at the end of the cycle.  Is this really needed? - YP
# @constraint(m, PM25[T+1] == PM25[1])
# @constraint(m, CO2[T+1] == CO2[1])
# @constraint(m, HUMID[T+1] == HUMID[1])

# maximum allowable indoor CO2 concentration [ppm].  Constraints do not apply to first timestep, since conditioning has not yet been applied.
@constraint(m, [t=2:T], CO2[t] <= CO2_MAX)

# maximum allowable indoor PM2.5 concentration [ug/m3].  Constraints do not apply to first timestep, since conditioning has not yet been applied.
@constraint(m, [t=2:T], PM25[t] <= PM25_MAX)

# maximum allowable indoor humidity ratio [gm water / gm of Dry Air].  Constraints do not apply to first timestep, since conditioning has not yet been applied.
@constraint(m, [t=2:T], HUMID[t] <= HUMID_MAX)

######################################
########### Print and solve ##########
######################################

#print(m)
solve(m)

CO2result = getvalue(CO2)
PM25result = getvalue(PM25)
HUMIDresult = getvalue(HUMID)

CMHresult = CMH
Nfilter = getvalue(N)
EnergyCost = getvalue(N)
PM25Absorbedresult = getvalue(roomPM25Absorption)
HUMIDAbsorbedresult = getvalue(kgMoistureRemoved)

fmax = maximum(CMHresult)
pmax = maximum(PM25Absorbedresult)
hmax = maximum(HUMIDAbsorbedresult)

cost = getobjectivevalue(m)

println("Total Cost [rmb]: ", cost)

println("Total Dehumidification Costs [RMB]: ", sum(HUMIDAbsorbedresult) * ISMRE * Celec)
println("Total Filter Costs [RMB]: ", Nfilter * Cfilter)
println("Total Fan Power Costs [RMB]: ", sum(CMH) * 365 * P * Celec)

println("Maximum Fan Load [cmh]: ", fmax)
println("Maximum PM2.5 Load [ug]: ", pmax)
println("Maximum Dehumidification Load [kg water]: ", hmax)

println("Nfilter [num]: ", Nfilter)


# Plotting -- use Pkg.add("Gadfly") if you don't have it installed
#= Each plot will work individually by plotting into the "plot" sector of Atom.
* To see them all together, make sure all the plots are defined and then run the
* vstack command in the REPL. That will make it open in your browser!
=#

CMH_plot =
    plot(
        x = 1:length(CMHresult),
        y = CMHresult,
        Geom.line,
        Guide.Title("CMH at timestep x")
    )

CO2_plot =
    plot(
        x = 0:length(CO2result) - 1,
        y = CO2result,
        Geom.line,
        Guide.Title("CO2 Concentration inside room (ppm)")
    )

PM25_plot =
    plot(
        x = 0:length(PM25result) - 1,
        y = PM25result,
        Geom.line,
        Guide.Title("PM2.5 Concentration (µg)")
    )

PM25_absorbed_plot =
    plot(
        x = 1:length(PM25Absorbedresult),
        y = PM25Absorbedresult,
        Geom.line,
        Guide.Title("PM2.5 Removed by Filters")
    )

HUMID_plot =
    plot(
        x = 0:length(HUMIDresult) - 1,
        y = HUMIDresult,
        Geom.line,
        Guide.Title("Humidity Ratio (kg Water / kg dry air)")
    )

HUMID_absorbed_plot =
    plot(
        x = 1:length(HUMIDAbsorbedresult),
        y = HUMIDAbsorbedresult,
        Geom.line,
        Guide.Title("Moisture Removed by Dehumidification")
    )


final = vstack(hstack(CMH_plot, CO2_plot), hstack(PM25_plot, PM25_absorbed_plot), hstack(HUMID_plot, HUMID_absorbed_plot))

img = SVG("Debugging Plots.svg", 12inch, 12inch)
draw(img, final)
