# Compare_JM_and_LM
Code from paper "Li, Wenhao, Liang Li, and Brad C. Astor. A comparison of two approaches to dynamic prediction: Joint modeling and landmark modeling. Statistics in Medicine (2023)".

'LM_Data_Generation.R' is the algorithm to generate data from landmark model.

'LM_fitting' is an example used to generate data from landmark model and using coxph to fit the data.
'example_data.RData' is an example data generated from landmark model. It can used to do dynamic prediction using file 'LM_DP_example'.

Standard package like 'JM' and 'joineRML' can be used for JM method including model fitting and prediction.
