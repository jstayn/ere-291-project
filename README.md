# ere-291-project
Group Project for ERE 291 in Julia for optimization of HVAC systems in a building in Shanghai.

NOTE: The rest of this is out of date as of 2018-03-19. See readme.txt in final_Schreiber_Stayner_Wang for a better readme as well as a cleaner version of this code. 



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


# Outer Loop Problem Statement:

Definitions: 
* Decision Variables
    * x[i] = Which FAU (Fresh Air Unit) you are selecting (binary) (x[2] == 1 implies FAU 2 is chosen)
    * y[j] = Whether a damper is installed in room j (binary)
* Parameters
    * c[i] = Cost of FAU_i [rmb]
    * c_d  = Cost of one damper [rmb]

Min f(x) = sum over i of (x[i] * c[i])  + sum over j of y[j]

Subject to:

* sum over i of (y[i]) == 1 (can only have one FAU)
* Possibly 0 < F (feasibility from operations regime)
* For all i, x[i] >= 0 and y[i] >= 0


# Inner Loop Problem Statement:

Definitions:

* Decision Variables
    * x[i] = Which FAU is selected (binary)
    * f[t] = fan output [m^3 / hr]
* Dependent Variables
    * N = number of filters required to achieve atmospheric constraints
    * n[t] = indoors concentration of CO2 at time t [ppm]
    * w[t] = indoors concentration of PM2.5 at time t [µg / m^3]
* Parameters
    * C_electricity = cost of electricity [rmb/kWh]
    * P[i] = relates CFM output to power input for AHU i (continuous parameter)
    * C_filter = cost per filter [rmb/filter]
    * W[t] = ambient concentration of PM2.5 at time t [µg / m^3]
    * S[t] = Source strength of CO2 at time t [cm^3 / hr]
    * R = efficiency of filter [% incoming PM2.5 absorbed]
    * N_max = maximum allowable concentration of CO2 [ppm]
    * W_max = maximum allowable concentration of PM 2.5 [µg / m^3]
    * k = maximum amount PM2.5 a filter can catch before being replaced [µg]

Min: C_var(x,f) = C_electricity * sum over i of (x[i] * P[i]) * sum over t of (f[t]) + C_filter * N

Subject to:

* For all t in [1, t_max]: n[t] == n[t-1] * S[t] - f[t] * (n[t-1] - n[t])
* For all t in [1, t_max]: w[t] == w[t-1] + f[t] * (R * W[t] - w[t-1])
* For all t in [1, t_max]: n[t] <= N_max
* For all t in [1, t_max]: w[t] <= N_max
* For all t in [1, t_max]: N = (sum over t of (f[t] * R * W[t])) / k 


