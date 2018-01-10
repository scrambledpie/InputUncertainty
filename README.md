## InputUncertainty
Code for the results produced in the paper "Bayesian Simulation Optimisation with Input Uncertainty, Michael Pearce, Juergen Branke. Winter Simulation Conference 2017".

# Upload_IU.R
R source code for the expriments, generate a test function with two inputs, collect evaluations from it in order to find the value of the oinput that maximises the average of the funtion over the second input.

# TotalData733
R results save file, open file with "Results = readRDS('TotalData733')", this contains 100 repetitions of each data collection method where each repetition has a different test function and noise value.

# Plot_Results
R source code to open and plot the TotalData733 file and reproduce the plots in the paper.
