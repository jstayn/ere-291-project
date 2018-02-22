# ere-291-project
Group Project for ERE 291 in Julia for optimization of HVAC systems in a building in Shanghai


# Notes 2018-02-22

Inner Loop/Outer Loop implementation?

Inner Loop: Best way to run certain machinery to produce optimum comfort with least cost

Outer Loop: Best machinery to buy for optimum comfort at the lowest cost. 

*Set desired comfort level and minimize cost* 

Outer loop inputs:
* Cost Data
* Feedback from Inner Loop
    * Feasibility
    * Energy usage to meet desired comfort levels


Outer Loop Problem statement:

Definitions: 
* Decision Variables
    * y[i] = Which Air Handler Unit you are selecting (binary) (y[2] == 1 implies AHU 2 is operating)
    * c[i] = Cost to run AHU i 
    * z[j] = Whether you want to buy dampers for room j (binary)
* Parameters
    * P[i] = relates CFM output to power input for AHU i (continuous parameter)

Max f(x) = sum over i of (y[i] * c[i]) * sum over i and j of (z[j] * C_dampers) + sum over t of (CFM_t * P[i] * y[i] * C_)

Subject to:

* sum over i of (y[i]) == 1 (for small room, only one AHU is on at a time)

