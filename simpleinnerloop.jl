PATH_TO_SOLVERS = ENV["ERE291_SOLVERS"]

using JuMP
using Clp
using Cbc
using AmplNLWriter

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
Cfilter = 5 #$/filter

k = 100000 #ug/filter
R = 0.8 #Efficiency of filter
T = 24

PM25_MAX = 100
CO2_MAX = 800
PM25_0 = 50
CO2_0 = AMB_CO2[1]

m = Model(solver=AmplNLSolver(joinpath(PATH_TO_SOLVERS,"knitro"), ["outlev=2"]))

@variable(m, CO2[1:T+1] >= 0)
@variable(m, PM25[1:T+1] >= 0)
@variable(m, CMH[1:T] >= 0)
@variable(m, N >= 0)

@objective(m, Min, Celec*P*sum(CMH) + Cfilter*N)

@constraint(m, [t=2:T+1], CO2[t] == CO2[t-1] + roomCO2Source[101][t-1]/roomData[101][2] - CMH[t-1]/(roomData[101][2])*(CO2[t-1]-AMB_CO2[t-1]))
@constraint(m, [t=2:T+1], PM25[t] == PM25[t-1] + CMH[t-1]/(roomData[101][2])*(R*AMB_PM2_5[t-1] - PM25[t-1]))
@constraint(m, PM25[1] == PM25_0)
@constraint(m, CO2[1] == CO2_0)
@constraint(m, PM25[T+1] == PM25[1])
@constraint(m, CO2[T+1] == CO2[1])
@constraint(m, [t=1:T+1], PM25[t] <= PM25_MAX)
@constraint(m, [t=1:T+1], CO2[t] <= CO2_MAX)
@NLexpression(m, N == R/k*sum(CMH[t]*PM25[2:end][t] for t in 1:T))

print(m)
solve(m)

CO2result = getvalue(CO2)
PM25result = getvalue(PM25)
CMHresult = getvalue(CMH)
Nfilter = getvalue(N)

fmax = maximum(CMHresult)
cost = getobjectivevalue(m)
