# InputUncertainty
Code for the results produced in the paper "Bayesian Simulation Optimisation with Input Uncertainty, Michael Pearce, Juergen Branke. Winter Simulation Conference 2017".

## BayesOpt_InputUncertainty.R
R source code for the expriments, generate a test function with two inputs, collect evaluations from it in order to find the value of the oinput that maximises the average of the funtion over the second input.

## TotalData733
R results save file, open file with "Results = readRDS('TotalData733')", this contains 100 repetitions of each data collection method where each repetition has a different test function and noise value. Each of 5 methods was applied 100 times each using two different input uncertainty distributions, uniform, and triangular. The file is a list with 1000 elements, 2 distributions * 5 algoroithms * 100 repetitions = 1000 experiments.

## Plot_Results.R
R source code to open and plot the TotalData733 file and reproduce the plots in the paper.

