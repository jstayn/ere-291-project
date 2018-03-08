# ERE 291 Project
# Simple Inner Loop
# Authors: Sam Schreiber, John Stayner, Yan-Ping Wang


####################################
######### Initialize tools #########
####################################

PATH_TO_SOLVERS = ENV["ERE291_SOLVERS"]

using JuMP
using Clp
using Cbc
using AmplNLWriter

using Gadfly



###############################
######### Define Sets #########
###############################

# Currently testing a period of 24 hours
T = 24

###################################################
############ Define parameters and data ###########
###################################################

data = readcsv("Master-Data_v2.0.csv")
AQdata = readcsv("AirQualityData20170301.csv")
rooms = data[2:end,3]

massPerson = 150 #lbs
met = 0.843318 #kcal/hr-lb
qCO2 = 241 #cm3/kcal

CO2pp = massPerson * met * qCO2 #cm3 CO2 / person-hr

roomData = Dict{Int64,Array{Float64}}()
roomOccupancy = Dict{Int64, Array{Int64}}() #Number of People in each hour
roomCO2Source = Dict{Int64, Array{Float64}}() #cm3 of CO2 Produced in each hour
for i = 1:length(rooms)
    roomData[rooms[i]] = data[i+1,4:6]
    roomOccupancy[rooms[i]] = data[i+1,7:30]
    roomCO2Source[rooms[i]] = CO2pp .* roomOccupancy[rooms[i]]
end

AMB_PM2_5 = AQdata[3:end,2] #ug/m3
AMB_CO2 = AQdata[3:end,3] #ppm

Celec = .2 #$/kWh
P = .2 #kWh/m3
Cfilter = 50 #$/filter

k = 100000 #ug/filter
R = 0.8 #Efficiency of filter

PM25_MAX = 35
CO2_MAX = 800
PM25_0 = 0
CO2_0 = AMB_CO2[1]


####################################
######### Initialize Model #########
####################################

m = Model(solver=AmplNLSolver(joinpath(PATH_TO_SOLVERS,"knitro"), ["outlev=1", "ms_enable=1"]))

####################################
######## Decision variables ########
####################################

@variable(m, CO2[1:T+1] >= 0, start=400)
@variable(m, PM25[1:T+1] >= 0, start=100)
@variable(m, CMH[1:T] >= 0, start=1)
# @variable(m, N >= 0, start=1)
@variable(m, roomPM25Sink[1:T] >= 0)
@NLexpression(m, N, R/k*sum(roomPM25Sink[t] for t in 1:T))

######################################
######## Objective Functions #########
######################################

# Minimize the total cost, equal to the sum of the cost of fan electricity over the operational period, plus the cost of PM2.5 filter replacements
@NLobjective(m, Min, Celec*P*sum(CMH[t] for t in 1:T) + Cfilter*N)

######################################
############# Constraints ############
######################################

# CO2 concentration in the room at time t is equal to the CO2 concentration in the room at t-1 plus (the CO2 mass introduced at t minus the CO2 mass removed at t) divided by the room volume
# The constraint is structured this way because the HVAC system can only react to rising CO2 levels at time t, it cannot pre-emptively condition the air at t-1
@NLconstraint(m, [t=2:T], CO2[t] == CO2[t-1] + (roomCO2Source[101][t] - CMH[t]*(CO2[t] - AMB_CO2[t]))/roomData[101][2])

# PM2.5 concentration in the room at time t is equal to the PM 2.5 in the room at t-1 plus the mass of PM2.5 introduced at time t minus the mass filtered at time t divided by the room volume
@constraint(m, [t=2:T], PM25[t] == PM25[t-1] + (CMH[t]*AMB_PM2_5[t] - CMH[t]*PM25[t] - roomPM25Sink[t])/roomData[101][2])

# set dummy initial condition to ambient concentrations
@constraint(m, PM25[1] == PM25_0)
@constraint(m, CO2[1] == CO2_0)

# CO2 and PM2.5 concentrations reset at the end of the cycle.  Is this really needed? - YP
@constraint(m, PM25[T+1] == PM25[1])
@constraint(m, CO2[T+1] == CO2[1])

# CO2 and PM2.5 levels must not exceed maximum allowable levels.  Constraints do not apply to first timestep, since conditioning has not yet been applied.
@constraint(m, [t=2:T], PM25[t] <= PM25_MAX)
@constraint(m, [t=2:T], CO2[t] <= CO2_MAX)


######################################
########### Print and solve ##########
######################################

#print(m)
solve(m)

CO2result = getvalue(CO2)
PM25result = getvalue(PM25)
CMHresult = getvalue(CMH)
Nfilter = getvalue(N)
PM25Absorbedresult = getvalue(roomPM25Sink)

fmax = maximum(CMHresult)
cost = getobjectivevalue(m)

println("Total Cost [rmb]: ", cost)
println("CMHresult [cmh]: ", CMHresult)
println("Maximum Capacity [cmh]: ", fmax)
println("CO2result [ppm]: ", CO2result)
println("PM25result [ug/m3]: ", PM25result)
println("PM25absorbedresult [ug]: ", PM25Absorbedresult)
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

final = vstack(hstack(CMH_plot, CO2_plot), hstack(PM25_plot, PM25_absorbed_plot))

img = SVG("Debugging Plots.svg", 12inch, 12inch)
draw(img, final)
