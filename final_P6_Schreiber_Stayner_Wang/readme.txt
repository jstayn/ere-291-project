ERE 291 Final Code Submission
Sam Schreiber, John Stayner, and Yan-Ping Wang

Readme.txt

The five top-level .jl files here are:
* LP_outer_loop.jl -- runs the Stage 2 outer loop
    * Dependencies: module_NLP_loop.jl, module_innerloop.jl
    * Can run 1 year in a reasonable amount of time (an hour or two)
* NLP_inner_loop.jl -- runs the Stage 1 inner loop. Note that we haven't implemented an outer loop for specifically for this one, but LP_outer_loop could easily be modified to do so.
    * Can run 1 week in a reasonable amount of time (a few hours)
    * Have not waited for 1 year to run yet.
* plot_LP_results.jl -- plots the results of running Stage 2 design on optimal AHU for one week
* plot_NLP_results.jl -- plots the results of running Stage 1 design on optimal AHU for one week
* NLP-sensitivity.jl -- runs a sensitivity analysis on starting points for Stage 1 design, outputs data to "sensitivity_output.csv"

There is one top-level .m file (Matlab):
* sensitivity_plot.m -- uses "sensivity_output.csv" to make a 3D surface plot of
                        sensitivity

There is one top-level .Rmd file (R Markdown):
* plot_pareto_curves.Rmd -- takes in a csv "pareto_data.csv" which is output by
                            LP_outer_loop.jl and generates graphs on Capital
                            Cost vs NPC as well as a simple sensitivity analysis
                            on discount rate.

The rest of the files are .csv inputs or outputs to the programs.

Thanks for a great quarter!
