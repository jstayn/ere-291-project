####################################
######### Initialize tools #########
####################################

PATH_TO_SOLVERS = ENV["ERE291_SOLVERS"]

using JuMP
using Clp
using Cbc
using AmplNLWriter

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
PM25_0 = 100
CO2_0 = AMB_CO2[1]
####################################
######### Initialize Model #########
####################################

m = Model(solver=AmplNLSolver(joinpath(PATH_TO_SOLVERS,"knitro"), ["outlev=2"]))

####################################
######## Decision variables ########
####################################


@variable(m, CO2[1:T+1] >= 0, start=400)
@variable(m, PM25[1:T+1] >= 0, start=100)
@variable(m, CMH[1:T] >= 0, start=1)
@variable(m, N >= 0, start=1)
@variable(m, A[1:T] >= 0)

######################################
######## Objective Functions #########
######################################

#Minimize the total cost, equal to the sum of the cost of fan electricity over the operational period, plus the cost of PM2.5 filter replacements
@objective(m, Min, Celec*P*sum(CMH[t] for t in 1:T) + Cfilter*N)

######################################
############# Constraints ############
######################################

#CO2 concentration in the room at time t is equal to the CO2 concentration in the room at t-1 plus (the CO2 mass introduced at t-1 minus the CO2 mass removed at t-1) divided by the room volume
#Consider concentrations to be measured at the END of each hour
@NLconstraint(m, [t=2:T], CO2[t] == CO2[t-1] + (roomCO2Source[101][t] - CMH[t]*(CO2[t]-AMB_CO2[t]))/roomData[101][2])

#PM2.5 concentration in the room at time t is equal to the PM 2.5 in the room at t-1 plus the mass of PM2.5 introduced at t-1 divided by the room volume
@constraint(m, [t=2:T], PM25[t] == PM25[t-1] + (CMH[t]*AMB_PM2_5[t] - A[t])/roomData[101][2])

#set dummy initial condition to ambient concentrations
@constraint(m, PM25[1] == PM25_0)
@constraint(m, CO2[1] == CO2_0)

#CO2 and PM2.5 concentrations reset at the end of the cycle.  Is this really needed? - YP
@constraint(m, PM25[T+1] == PM25[1])
@constraint(m, CO2[T+1] == CO2[1])

#CO2 and PM2.5 levels must not exceed maximum allowable levels.  Constraints do not apply to first timestep, since conditioning has not yet been applied.
@constraint(m, [t=2:T], PM25[t] <= PM25_MAX)
@constraint(m, [t=2:T], CO2[t] <= CO2_MAX)

@constraint(m, N == R/k*sum(A[t] for t in 1:T))

######################################
########### Print and solve ##########
######################################

print(m)
solve(m)

CO2result = getvalue(CO2)
PM25result = getvalue(PM25)
CMHresult = getvalue(CMH)
Nfilter = getvalue(N)
PM25Absorbedresult = getvalue(A)

fmax = maximum(CMHresult)
cost = getobjectivevalue(m)

println("Total Cost [rmb]: ", cost)
println("CMHresult [cmh]: ", CMHresult)
println("Maximum Capacity [cmh]: ", fmax)
println("CO2result [ppm]: ", CO2result)
println("PM25result [ug/m3]: ", PM25result)
println("PM25absorbedresult [ug]: ", PM25Absorbedresult)
println("Nfilter [num]: ", Nfilter)
