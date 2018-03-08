# ERE 291 Project
# Simple Outer Loop (find best design constraints)
# Authors: Sam Schreiber, John Stayner, Yan-Ping Wang
# 2018-03-08


# Simple outer loop
using inner_loop


T = 24
P = 0.5
capacity = 10000

equip_data = readcsv("Equipment_Costs.csv")[2:end,:]

equip_costs = equip_data[:, 2]
equip_capacity = equip_data[:, 3]
equip_efficiency = equip_data[:, 4] / 100

feasibility, cost = operations_optimization(T, P, capacity)

optimal_cost = Inf
best = 0


for i = 1:size(equip_data)[1]
    feasibility, ops_cost = operations_optimization(T, 1 / equip_efficiency[i],
                                                    equip_capacity[i])
    if feasibility
        total_cost =  365 * ops_cost + equip_costs[i]
        if total_cost < optimal_cost
            optimal_cost = total_cost
            best = i
        end
    end
end

if best == 0
    best = "Infeasible"
end

println("Best AHU: ", best, "\n",
        "Optimal Cost: ", optimal_cost)
