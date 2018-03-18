# ERE 291 Project
# Simple Outer Loop (find best design constraints)
# Authors: Sam Schreiber, John Stayner, Yan-Ping Wang
# 2018-03-08


# Simple outer loop
using inner_loop
using Gadfly

T = 8760

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

capexResults = []
opexResults = []
NPCResults = []
feasibilityResults = []

fmaxResults = []
pmaxResults = []
hmaxResults = []

NLPstatusResults = []

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

for i = 1:size(equip_data)[1]

    for j = 1:length(damperPlacementOptions)

        feasibility, opex, fmax, pmax, hmax, NLPstatus = OpsOptLP(T, equip_efficiency[i],
                                                        equip_capacity[i], ISMRE[i], damperPlacementOptions[j])
        capex = equip_costs[i] + sum(damper_cost * (1 - damperPlacementOptions[j][k]) for k = 1:length(damperPlacementOptions[j]))
        NPC = capex + sum(opex/(1+discountRate)^p for p=1:nYears)

        push!(capexResults, capex)
        push!(opexResults, opex)
        push!(NPCResults, NPC)
        push!(feasibilityResults, feasibility)

        push!(fmaxResults, fmax)
        push!(pmaxResults, pmax)
        push!(hmaxResults, hmax)

        push!(NLPstatusResults, NLPstatus)

        println("AHU Number ", i, " and damper Option ", j, " Results:")
        println("Capex: ", capex)
        println("Opex: ", opex)
        println("NPC: ", NPC)
        println("Feasibility: ", feasibility)
    end
end

writecsv("capex.csv", capexResults)
writecsv("opex.csv", opexResults)
writecsv("NPC.csv", NPCResults)
writecsv("feasibility.csv", feasibilityResults)
writecsv("fmax.csv", fmaxResults)
writecsv("pmax.csv", pmaxResults)
writecsv("hmax.csv", hmaxResults)
#writecsv("NLPStatus.csv", NLPstatusResults)

println("NLPStatus: ", NLPstatusResults)

#ops_vs_cap = plot(
    #x = capexResults,
    #y = NPCResults,
    #Geom.point,
    #Guide.Title("Figure 3: Relationship between\nNPC and Capital Costs"),
    #Guide.XLabel("Capital Costs (RMB)"),
    #Guide.YLabel("NPC (RMB)")
    #)

#img = SVG("NPC vs Capital Costs.svg", 4inch, 4inch)
#draw(img, NPC_vs_cap)
