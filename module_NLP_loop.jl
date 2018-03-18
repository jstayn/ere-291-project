# Working CMH Solver with dampers to be called by module_innerloop.jl

module NLP_loop

using JuMP
using Clp
using Cbc
using AmplNLWriter
using Gadfly

export CMH_opt

####################################
######## Initialize packages #######
####################################

PATH_TO_SOLVERS = ENV["ERE291_SOLVERS"]

function CMH_opt(T, N, damperPlacement)

    data = readcsv("Room-Data_v3.0.csv")
    AQdata = readcsv("AirQualityData2016.csv")
    rooms = data[2:end,3]

    ###################################################
    ############ Define parameters and data ###########
    ###################################################

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
    roomVolume = Dict{Int64,Float64}()
    roomDiffusers = Dict{Int64,Int64}()
    roomOccupancy = Dict{Int64, Array{Int64}}() #number of People in each hour
    roomCO2Source = Dict{Int64, Array{Float64}}() #cm3 of CO2 Produced in each hour

    for i = 1:N
        roomData[rooms[i]] = data[i+1,4:7]
        roomVolume[rooms[i]] = data[i+1,5]
        roomDiffusers[rooms[i]] = data[i+1,7]
        roomOccupancy[rooms[i]] = data[i+1,8:31]
        roomCO2Source[rooms[i]] = CO2pp .* roomOccupancy[rooms[i]]
    end

    # Number of air diffusers per room
    diffusers = data[2:end, 7]
    numDiffusers = sum(diffusers)

    # initial indoor CO2 concentration at t = 0
    CO2_0 = AMB_CO2[1]

    # efficiency of the fan of FAU model i [kWh / m3 air]
    P = .00045

    ### Cost Parameters ###

    # cost of electricity at time t [RMB/kWh]
    Celec = 1.4


    ####################################
    ######### Initialize Model #########
    ####################################

    m = Model(solver=AmplNLSolver(joinpath(PATH_TO_SOLVERS,"knitro"), ["outlev=0"])) #, "ms_enable=1"]))


    ####################################
    ######## Decision variables ########
    ####################################

    # rate of air intake at time t [m3 / hr]
    @variable(m, CMHCentral[1:T] >= 0, start=1)

    # CO2 concentration in room i at time t [m3 / ug]
    @variable(m, CO2[1:N, 1:T] >= 0, start=1)

    # Damper positions in room i at time t [%]
    @variable(m, 0 <= damperPosition[1:N, 1:T] <= 1, start=1)

    # Damper positions in room i at time t [%]
    @variable(m, CMHRoom[1:N, 1:T] >= 0, start=1)


    ######################################
    ######## Dependent Variables #########
    ######################################

    # force the damper open if damper_placement == 0
    # If damper exists (damperPlacement[i] == 0), damperPosition is unaffected.
    # If damper does not exist (damperPlacement[i] == 1), then damperPosition == 1.
    @NLexpression(m, damperPositionFinal[i=1:N, t=1:T], min(damperPosition[i,t] + damperPlacement[i], 1))


    ######################################
    ######## Objective Functions #########
    ######################################

    # Minimize the total cost, equal to the sum of the cost of fan electricity over the operational period, plus the cost of PM2.5 filter replacements
    @objective(m, Min, P * Celec * sum(CMHCentral[t] for t=1:T))


    ######################################
    ############# Constraints ############
    ######################################

    # CO2 concentration in the room at time t is equal to the CO2 concentration in the room at t-1 plus (the CO2 mass introduced at t minus the CO2 mass removed at t) divided by the room volume
    # The constraint is structured this way because the HVAC system can only react to CO2 levels at time t, it cannot pre-emptively condition the air at t-1
    @NLconstraint(m, [i=1:N, t=2:T], CO2[i,t] >= CO2[i,t-1] + (roomCO2Source[rooms[i]][t] + CMHRoom[i,t]*AMB_CO2[t] - CMHRoom[i,t]*CO2[i,t])/roomVolume[rooms[i]])

    # set  initial condition to ambient concentrations
    @constraint(m, [i=1:N], CO2[i,1] == CO2_0)

    # maximum allowable indoor CO2 concentration [ppm].  Constraints do not apply to first timestep, since conditioning has not yet been applied.
    @constraint(m, [i=1:N, t=2:T], CO2[i,t] <= CO2_MAX)

    # the CMH in each room is equal to the total CMH divided by the ratio of the damper opening to the current room over the total damper openings.  Assumes all ducts are the same size
    #@NLconstraint(m, [i=1:N, t=1:T], CMHRoom[i,t] <= CMHCentral[t] * (damperPosition[i,t] / sum(damperPosition[i,t] for i=1:N)))
    @NLconstraint(m, [i=1:N, t=1:T], CMHRoom[i,t] <= CMHCentral[t] * (damperPositionFinal[i,t]*diffusers[i]) / (sum(damperPositionFinal[i,t] * diffusers[i] for i in 1:N)))

    ######################################
    ########### Print and solve ##########
    ######################################

    #print(m)
    solve(m)

    CO2result = getvalue(CO2)
    CO2resultSum = sum(CO2result, 1)

    CMHresult = getvalue(CMHCentral)
    CMHRoomResult = getvalue(CMHRoom)

    fmax = maximum(CMHresult)

    cost = getobjectivevalue(m)

    println("Total Cost [rmb]: ", cost)
    println("Total Fan Power Costs [RMB]: ", sum(CMHresult) * P * Celec)
    println("Maximum Fan Load [cmh]: ", fmax)

    CMHRoomResultRounded = round.(Int, CMHresult)

    return CMHRoomResult

end

end
