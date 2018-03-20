####################################
######## Initialize packages #######
####################################

PATH_TO_SOLVERS = ENV["ERE291_SOLVERS"]

using JuMP
using AmplNLWriter

function innerloop(x, y, T=24, P=.00045)

    ###################################################
    ############ Define parameters and data ###########
    ###################################################

    ### Air Quality Parameters ###

    data = readcsv("Room-Data_v3.0.csv")
    AQdata = readcsv("AirQualityData2016.csv")
    rooms = data[2:end,3]
    N = length(rooms)

    # Ambient Air Quality Parameters
    AMB_CO2 = AQdata[3:end,3] #[ppm]
    AMB_PM2_5 = AQdata[3:end,2] #[ug/m3]
    AMB_HUMID = AQdata[3:end,4] #[gm water / gm dry air]

    # Room Air Quality Maximum Tolerances
    CO2_MAX = 1000 #[ppm]
    PM25_MAX = 35 #[ug/m3]
    HUMID_MAX = 0.012 #[gm water / gm dry air]

    # Room CO2 source per person
    massPerson = 150 #lbs
    met = 0.843318 #kcal/hr-lb
    qCO2 = 241 #cm3/kcal
    CO2pp = massPerson * met * qCO2 #cm3 CO2 / person-hr

    # Room humidity source per person
    humidityPerPerson = 0.15 #[kg / person-hr]

    # Room data dicts
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

    # Number of air diffusers per room
    diffusers = data[2:end, 7]
    numDiffusers = sum(diffusers)

    # pressure drop through system
    pressureDropDucts = 103 #[Pa]
    pressureDropHEPAFilter = 250 #[Pa]
    pressureDropSystem = pressureDropDucts + pressureDropHEPAFilter

    # initial indoor conditions
    CO2_0 = AMB_CO2[1]
    PM25_0 = AMB_PM2_5[1]
    HUMID_0 = AMB_HUMID[1]

    # density of air [kg / m3]
    airDensity = 1.225

    ### Equipment Parameters ###

    # efficiency of FAU fan
    P = .00045 #[kWh / m3 air]

    # PM2.5 absorption capacity per filter before replacement
    k = 10^7 #[ug/filter]

    # minimum integrated seasonal moisture removal efficiency (ISMRE)
    ISMRE = 0.5 #[kWh/kg]

    # efficiency of PM2.5 filter
    R = 0.8

    ### Cost Parameters ###

    # cost of electricity at time t [RMB/kWh]
    Celec = 1.4

    # cost per in-room PM2.5 filter [RMB/filter]
    Cfilter = 50

    ####################################
    ######### Initialize Model #########
    ####################################


    #m = Model(solver=AmplNLSolver(joinpath(PATH_TO_SOLVERS,"knitro"), ["outlev=0", "ms_enable=1"]))
    m = Model(solver = AmplNLSolver(joinpath(PATH_TO_SOLVERS,"knitro"), ["outlev=2", "maxit=200"]))

    ####################################
    ######## Decision variables ########
    ####################################

    # rate of PM2.5 absorption in room i at time t [ug / hr]
    @variable(m, roomPM25Absorption[1:N, 1:T] >= 0)

    # rate of humidity removal at time t [kg / hr]
    @variable(m, kgMoistureRemoved[1:N, 1:T] >= 0)

    # CO2 concentration in room i at time t [m3 / ug]
    @variable(m, CO2[1:N, 1:T] >= 0, start=y)

    # PM2.5 concentration in room i at time t [ppm]
    @variable(m, PM25[1:N, 1:T] >= 0)

    # humidity ratio of indoor air at time t [gm water / gm of dry air]
    @variable(m, HUMID[1:N, 1:T] >= 0)

    # rate of air intake at time t [m3 / hr]
    @variable(m, CMHCentral[1:T] >= 0, start=x)

    # damper positions in room i at time t [%]
    @variable(m, 0 <= damperPosition[1:N, 1:T] <= 1)

    @variable(m, CMHRoom[1:N, 1:T] >= 0)

    #####################################
    ######## Dependent variables ########
    #####################################

    # number of PM2.5 filters required for all rooms
    @expression(m, numFilters, sum(roomPM25Absorption[i,t] for i in 1:N for t in 1:T) / k)

    # rate of incoming humidity from FAU at time t [kg / hr]
    @NLexpression(m, kgMoistureFAUIn[i=1:N, t=1:T], CMHRoom[i,t] * airDensity * AMB_HUMID[t])

    # rate of outgoing humidity from FAU at time t [kg / hr]
    @NLexpression(m, kgMoistureFAUOut[i=1:N, t=1:T], CMHRoom[i,t] * airDensity * HUMID[i,t])

    ######################################
    ######## Objective Functions #########
    ######################################

    # Minimize the total cost, equal to the sum of the cost of fan electricity over the operational period, plus the cost of PM2.5 filter replacements
    @objective(m, Min, Celec*((P + pressureDropSystem / 3600/ 1000) * sum(CMHCentral[t] for t in 1:T) + ISMRE*sum(kgMoistureRemoved[i,t] for i in 1:N for t in 1:T)) + Cfilter*numFilters)

    ######################################
    ############# Constraints ############
    ######################################

    @NLconstraint(m, [i=1:N, t=1:T], CMHRoom[i,t] == CMHCentral[t] * (damperPosition[i,t]*diffusers[i]) / (sum(damperPosition[j,t] * diffusers[j] for j in 1:N)))

    # CO2 concentration in the room at time t is equal to the CO2 concentration in the room at t-1 plus (the CO2 mass introduced at t minus the CO2 mass removed at t) divided by the room volume
    # The constraint is structured this way because the HVAC system can only react to CO2 levels at time t, it cannot pre-emptively condition the air at t-1
    # What if we just assumed an air change rate based on the maximum CO2
    @NLconstraint(m, [i=1:N, t=2:T], CO2[i,t] >= CO2[i,t-1] + (roomCO2Source[rooms[i]][(t-1)%24 + 1] + CMHRoom[i,t]*(AMB_CO2[t] - CO2[i,t]))/roomData[rooms[i]][2])

    # PM2.5 concentration in the room at time t is equal to the PM 2.5 in the room at t-1 plus the mass of PM2.5 introduced at time t minus the mass filtered at time t divided by the room volume
    @NLconstraint(m, [i=1:N, t=2:T], PM25[i,t] >= PM25[i,t-1] + (CMHRoom[i,t]*(AMB_PM2_5[t] - PM25[i,t]) - roomPM25Absorption[i,t])/roomData[rooms[i]][2])

    # Humidity ratio in the room at time t is equal to the Humidity ratio  in the room at t-1 plus (the H2O mass introduced at t minus the H2O mass removed at t) divided by the room volume
    # The constraint is structured this way because the HVAC system can only react to humidity levels at time t, it cannot pre-emptively condition the air at t-1
    @NLconstraint(m, [i=1:N, t=2:T], HUMID[i,t] >= HUMID[i,t-1] + (roomHumidSource[rooms[i]][(t-1)%24 + 1] + kgMoistureFAUIn[i,t] - kgMoistureFAUOut[i,t] - kgMoistureRemoved[i,t])/(roomData[rooms[i]][2] * airDensity))

    # set  initial condition to ambient concentrations
    @NLconstraint(m, [i=1:N], PM25[i,1] == PM25_0)
    @NLconstraint(m, [i=1:N], CO2[i,1] == CO2_0)
    @NLconstraint(m, [i=1:N], HUMID[i,1] == HUMID_0)

    # CO2, PM2.5 and Humidity concentrations reset at the end of the cycle.
    # @NLconstraint(m, [i=1:N], PM25[i,T] == PM25[i,1])
    # @NLconstraint(m, [i=1:N], CO2[i,T] == CO2[i,1])
    # @NLconstraint(m, [i=1:N], HUMID[i,T] == HUMID[i,1])

    # maximum allowable indoor CO2 concentration [ppm].  Constraints do not apply to first timestep, since conditioning has not yet been applied.
    @NLconstraint(m, [i=1:N, t=2:T], CO2[i,t] <= CO2_MAX)

    # maximum allowable indoor PM2.5 concentration [ug/m3].  Constraints do not apply to first timestep, since conditioning has not yet been applied.
    @constraint(m, [i=1:N, t=2:T], PM25[i,t] <= PM25_MAX)

    # maximum allowable indoor humidity ratio [gm water / gm of Dry Air].  Constraints do not apply to first timestep, since conditioning has not yet been applied.
    @constraint(m, [i=1:N, t=2:T], HUMID[i,t] <= HUMID_MAX)

    ######################################
    ########### Print and solve ##########
    ######################################

    status = solve(m)

    if status == :Optimal res = 1
    elseif status == :Unbounded res = 2
    elseif status == :Infeasible res = 3
    elseif status == :UserLimit res = 4
    elseif status == :Error res = 5
    elseif status == :NotSolved res = 6
    end

    CHMroomoutput = getvalue(CMHRoom)
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

    Nfilter = getvalue(numFilters)


    PM25Absorbedresult = sum(getvalue(roomPM25Absorption),1)
    HUMIDAbsorbedresult = sum(getvalue(kgMoistureRemoved),1)
    CMHresult = getvalue(CMHCentral)

    fmax = maximum(CMHresult)
    pmax = maximum(PM25Absorbedresult)
    hmax = maximum(HUMIDAbsorbedresult)

    cost = getobjectivevalue(m)

    println("Total Cost [rmb]: ", cost)
    println("Total Dehumidification Costs [RMB]: ", sum(HUMIDAbsorbedresult) * ISMRE * Celec)
    println("Total Filter Costs [RMB]: ", Nfilter * Cfilter)
    println("Total Fan Power Costs [RMB]: ", sum(CMHresult) * P * Celec)
    println("Maximum Fan Load [cmh]: ", fmax)
    println("Maximum PM2.5 Load [ug]: ", pmax)
    println("Maximum Dehumidification Load [kg water]: ", hmax)
    println("Nfilter [num]: ", Nfilter)

    return [x y cost res]
end

CMHmin = 0
CMHmax = 10000
CO2min = 0
CO2max = 1000
N = 21

CMHstep = (CMHmax-CMHmin)/(N-1)
CO2step = (CO2max-CO2min)/(N-1)

results = zeros(N*N,4)
k=1

for i = 1:N
    for j = 1:N
        results[k,:] = innerloop((i-1)*CMHstep,(j-1)*CO2step)
        k += 1
    end
end

 filename = "sensitivity-output.csv"
 writecsv(filename, results)
