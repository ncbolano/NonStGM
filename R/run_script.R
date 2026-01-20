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
JJ = variable_list[1] ; p = variable_list[2] ; nu = variable_list[3]

