# Test for ERE 291

module inner_loop

using JuMP
using Clp
using Cbc
using AmplNLWriter
using Gadfly

export operations_optimization
export one_day
####################################
######### Initialize tools #########
####################################

PATH_TO_SOLVERS = ENV["ERE291_SOLVERS"]

function operations_optimization(T, P, capacity)


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
    #P = .2 #kWh/m3
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

    m = Model(solver=AmplNLSolver(joinpath(PATH_TO_SOLVERS,"knitro"), ["outlev=0", "ms_enable=1"]))

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
    println("Solved with P = ", P, ", and f_max = ", f_max)

    return fmax <= capacity, cost

end

function one_day(P)

    T = 24 # Hours

    ###################################################
    ############ Define parameters and data ###########
    ###################################################

    ### Air Quality Parameters ###

    data = readcsv("Master-Data_v2.0.csv")
    AQdata = readcsv("AirQualityData2016.csv")
    rooms = data[2:end,3]
    N = length(rooms)

    # CO2 concentration of outdoor air at time t [ppm]
    AMB_CO2 = AQdata[3:end,3]

    # maximum allowable indoor CO2 concentration [ppm]
    CO2_MAX = 800

    # Calculate rate of CO2 emissions from occupancy at time t [cm3 / hr]
    massPerson = 150 #lbs
    met = 0.843318 #kcal/hr-lb
    qCO2 = 241 #cm3/kcal

    CO2pp = massPerson * met * qCO2 #cm3 CO2 / person-hr

    roomData = Dict{Int64,Array{Float64}}()
    roomOccupancy = Dict{Int64, Array{Int64}}() #number of People in each hour
    roomCO2Source = Dict{Int64, Array{Float64}}() #cm3 of CO2 Produced in each hour
    roomCO2Max = Dict{Int64, Float64}() #cm3 of CO2 Produced in each hour
    roomHumidSource = Dict{Int64, Array{Float64}}() #kg H20 produced per hour

    for i = 1:N
        roomData[rooms[i]] = data[i+1,4:6]
        roomOccupancy[rooms[i]] = data[i+1,7:30]
        roomCO2Source[rooms[i]] = CO2pp .* roomOccupancy[rooms[i]]
    end

    # initial indoor CO2 concentration at t = 0
    CO2_0 = CO2_MAX

    # density of air [kg / m3]
    airDensity = 1.225

    ### Cost Parameters ###

    # cost of electricity at time t [RMB/kWh]
    Celec = 1.4

    # cost per in-room PM2.5 filter [RMB/filter]
    Cfilter = 50


    ####################################
    ######### Initialize Model #########
    ####################################

    m = Model(solver=AmplNLSolver(joinpath(PATH_TO_SOLVERS,"knitro"), ["outlev=2"]))

    ####################################
    ######## Decision variables ########
    ####################################

    # CO2 concentration in room i at time t [m3 / ug]
    @variable(m, CO2[1:N, 1:T] >= 0)
    @variable(m, CMH[1:N, 1:T] >= 0)

    ######################################
    ######## Objective Functions #########
    ######################################

    # Minimize the total cost, equal to the sum of the cost of fan electricity over the operational period, plus the cost of PM2.5 filter replacements
    @NLobjective(m, Min, Celec*P*sum(CMH[i] for i = 1:N*T))


    ######################################
    ############# Constraints ############
    ######################################

    # CO2 concentration in the room at time t is equal to the CO2 concentration in the room at t-1 plus (the CO2 mass introduced at t minus the CO2 mass removed at t) divided by the room volume
    # The constraint is structured this way because the HVAC system can only react to CO2 levels at time t, it cannot pre-emptively condition the air at t-1
    # What if we just assumed an air change rate based on the macapacityimum CO2
    @NLconstraint(m, [i=1:N, t=2:T], CO2[i,t] >= CO2[i,t-1] + (roomCO2Source[rooms[i]][(t-1)%24 + 1] + CMH[i,(t-1)%24 + 1]*(AMB_CO2[t] - CO2[i,t]))/roomData[rooms[i]][2])

    # set  initial condition to ambient concentrations
    @constraint(m, [i=1:N], CO2[i,1] == CO2_0)

    # CO2, PM2.5 and Humidity concentrations reset at the end of the cycle.
    @constraint(m, [i=1:N], CO2[i,T] == CO2[i,1])

    # maximum allowable indoor CO2 concentration [ppm].  Constraints do not apply to first timestep, since conditioning has not yet been applied.
    @constraint(m, [i=1:N, t=2:T], CO2[i,t] <= CO2_MAX)

    ######################################
    ########### Print and solve ##########
    ######################################

    #print(m)
    solve(m)

    return getvalue(CMH)


end

#CMH = one_day(0.4)

#plot(x = 1:24, y = CMHresult, Geom.line)
#plot(x = 1:24, y = sum(values(roomCO2Source), 1), Geom.line)
end
