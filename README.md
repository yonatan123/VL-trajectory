
Code for: "Estimating Mean Viral Load Trajectory from Intermittent Longitudinal Data and Unknown Time Origins"

This folder includes Matlab files for numerical simulations

Type 'syntheticexperiment' on command window to run the simualtions. 

Files included: 
syntheticexperiment.m: main script
synthparamteres.mat: true parameters for the sample generation as described in Section 5.
initializeparams.m: random initialization for the EM algorithm
generatesamples.m: sample generation according to the model introduces in the paper
gamma_approximation.m: curve approximation using Gamma curve
emfit.m: EM-algorithm without constraints
emfitunimodal.m: EM-algorithm under unimodality constraint
emfitgamma.m: EM-algorithm under Gamma constraint


Data structure: 
The Ct-value data, denoted in the code by D, is a struct object and should include the following fields:

ct2:  n x 2 double matrix, each row representing a pair of Ct-values, with columns 1,2 consisting of 1st and 2nd Ct-value, respectively
delta2: n x 1 double matrix (of integers) representing the time difference between measurements

(optional)
ct3:  n x 3 double matrix, each row representing a triplet of Ct-values, with columns 1,2,3 consisting of 1st, 2nd and 3rd Ct-value, respectively
delta2: n x 2 double matrix (of integers) representing the time difference between each of the two consecutive measurements

similarly ct3, delta3, ct4, delta4, etc. 
