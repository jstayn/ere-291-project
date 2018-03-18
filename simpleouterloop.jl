# ERE 291 Project
# Simple Outer Loop (find best design constraints)
# Authors: Sam Schreiber, John Stayner, Yan-Ping Wang
# 2018-03-08


# Simple outer loop
using inner_loop
using Gadfly


T = 24
P = 0.5
capacity = 10000

equip_data = readcsv("Equipment_Costs.csv")[2:end,:]

equip_costs = equip_data[:, 2]
equip_capacity = equip_data[:, 3]
equip_efficiency = equip_data[:, 4] / 100

optimal_cost = Inf
best = 0

ops_costs = []

# Note that right now our inner loop is implemented such that the efficiency
# number it takes is in kWh/m^3 where more efficiency numbers are given in
# m^3/kWh, so we take the inverse of the efficiency when we pass it in on line
# 33.

for i = 1:size(equip_data)[1]
    feasibility, ops_cost = operations_optimization(T, 1 / equip_efficiency[i],
                                                    equip_capacity[i])
    if feasibility
        total_cost =  ops_cost + equip_costs[i]
        push!(ops_costs, ops_cost)
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


function cap_costs_scale(x)
    scale_factor = 1000
    "\$$(x / scale_factor)"
end

function NPV_scale(y)
    scale_factor = 100
    "\$$(y / scale_factor)"
end


light_theme = Theme(
    background_color = "white",
    panel_fill = "white",
    default_color = "blue"
    )

ops_vs_cap = plot(
    x = equip_costs,
    y = ops_costs,
    Geom.point,
    Guide.Title("Figure 3: Relationship between\nOperation Costs and Capital Costs"),
    Guide.XLabel("Capital Costs (Thousands of RMB)"),
    Guide.YLabel("Operating Costs (RMB)"),
    Scale.x_continuous(labels = cap_costs_scale),
    Scale.y_continuous(labels = NPV_scale),
    light_theme
    )

img = SVG("Op Costs vs Capital Costs.svg", 4inch, 4inch)
draw(img, ops_vs_cap)
