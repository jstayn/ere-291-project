# Test for ERE 291

module inner_loop

include("module_NLP_loop.jl")

using JuMP
using Clp
using Cbc
using AmplNLWriter
using Gadfly


using inner_loop.NLP_loop

export OpsOptLP

####################################
######### Initialize tools #########
####################################

PATH_TO_SOLVERS = ENV["ERE291_SOLVERS"]

function OpsOptLP(T, P, capacity, ISMRE, damperPlacement)

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

    CMH = zeros(N,T)

    for i in 1:N
        for j in 1:T_repeat
            CMH[i,j] = max(CMHOrig[i,j], CMHOptResult[i,j])
        end
    end

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
    @objective(m, Min, Celec*((P + pressureDropSystem / 3600 / 1000)*sum(CMH)*T/24 + ISMRE*sum(kgMoistureRemoved)) + Cfilter*numFilters)


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

    CO2result = sum(getvalue(CO2),1)
    PM25result = sum(getvalue(PM25),1)
    HUMIDresult = sum(getvalue(HUMID),1)

    Nfilter = getvalue(numFilters)


    PM25Absorbedresult = sum(getvalue(roomPM25Absorption),1)
    HUMIDAbsorbedresult = sum(getvalue(kgMoistureRemoved),1)

    pmax = maximum(PM25Absorbedresult)
    hmax = maximum(HUMIDAbsorbedresult)

    cost = getobjectivevalue(m)

    println("Total Cost [rmb]: ", cost)

    println("Total Dehumidification Costs [RMB]: ", sum(HUMIDAbsorbedresult) * ISMRE * Celec)
    println("Total Filter Costs [RMB]: ", Nfilter * Cfilter)
    println("Total Fan Power Costs [RMB]: ", sum(CMH)*T/24 * P * Celec)

    println("Maximum Fan Load [cmh]: ", fmax)
    println("Maximum PM2.5 Load [ug]: ", pmax)
    println("Maximum Dehumidification Load [kg water]: ", hmax)

    return fmax <= capacity, cost, fmax, pmax, hmax, NLPstatus

end

end
