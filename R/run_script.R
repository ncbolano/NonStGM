# Complete run-through (req. functions)
source('extractVariables.R')
source('KernelWeights.R')
# source('Combined_MK_Estimation.R')
source('DFTransform.R')
source('R_Hat_Creation.R')
source('BetaFunctions.R')
source('VarianceFunctions.R')

# Extract key variables from seq. x
variable_list = extractVariables(x)
J = variable_list[1] ; p = variable_list[2] ; nu = variable_list[3]
coefnum = (2 * p * nu) + p - 1

# Gathering frequencies (k) + local smoothing values (m_i) for each frequency
k = extractK(J,coefnum)

# End of chosen list of two k frequencies
# Obtaining M_list for each k frequency
M_list = local_M_selection(JJ, k, M_grid)
M_grid = unique(round(seq(max(coefnum , n^(1/5)), n^(1/2), length.out = 50)))


CV_scores = numeric(length(M_grid))
M_list = local_M_selection(JJ, k, M_grid)






