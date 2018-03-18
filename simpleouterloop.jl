# ERE 291 Project
# Simple Outer Loop (find best design constraints)
# Authors: Sam Schreiber, John Stayner, Yan-Ping Wang
# 2018-03-08


# Simple outer loop
using inner_loop
using Gadfly


T = 8760
capacity = 10000

equip_data = readcsv("Equipment_Costs.csv")[2:end,:]
damper_data = readcsv("Damper_Data.csv")[2:end,:]

equip_costs = equip_data[:, 2]
equip_capacity = equip_data[:, 3]
equip_efficiency = equip_data[:, 4]
ISMRE = equip_data[:, 5]

damper_placement = damper_data[:, 2]
damper_costs = damper_data[:, 3]

optimal_cost = Inf
best = 0

capexResults = []
NPVResults = []
nYears = 20
discountRate = 0.06

feasibilityResults = []

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

    for j = 1:13

        feasibility, ops_cost = OpsOptLP(T, equip_efficiency[i],
                                                        equip_capacity[i], ISMRE[i], damperPlacementOptions[j])
        capex = equip_costs[i] + sum(damper_costs * damperPlacementOptions[j][k] for k = 1:12)
        push!(capexResults, capex)

        NPV = capex + sum(ops_cost/(1+discountRate)^p for p=1:nYears)
        push!(NPVResults, NPV)

        push!(feasibilityResults, feasibility)

        if feasibility
            if total_cost < optimal_cost
                optimal_cost = total_cost
                best = i
            end
        end
    end
end

if best == 0
    best = "Infeasible"
end

println("Best AHU: ", best, "\n",
        "Optimal Cost: ", optimal_cost)

writecsv("capex.csv", capexResults)
writecsv("NPV.csv", NPVResults)
writecsv("feasibility.csv", feasibilityResults)

ops_vs_cap = plot(
    x = capex,
    y = NPVResults,
    Geom.point,
    Guide.Title("Figure 3: Relationship between\nNPV and Capital Costs"),
    Guide.XLabel("Capital Costs (RMB)"),
    Guide.YLabel("NPV (RMB)")
    )

img = SVG("NPV vs Capital Costs.svg", 4inch, 4inch)
draw(img, NPV_vs_cap)
