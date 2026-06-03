
library(foreach)
library(doParallel)
library(parallel)
# ============================================================
# PARALLEL SIMULATION RUNNER FOR p=3
# ============================================================
run_NonStGM_simulation_p3_parallel = function(R = 100, alpha = 0.05, burnin = 20, m = 4096,
                                              TV_size = 0.6, nu = 2, L = 1,
                                              Kernel = 'Kernel_Triangular', seed = 0,
                                              n_cores = NULL) {

  if (is.null(n_cores)) {
    n_cores = detectCores() - 1
  }

  cl = makeCluster(n_cores)
  registerDoParallel(cl)

  clusterEvalQ(cl, {
    source('extractVariables.R')
    source('exampleVAR.R')
    source('KernelWeights.R')
    source('Combined_MK_Estimation.R')
    source('DFTransform.R')
    source('R_Hat_Creation.R')
    source('BetaFunctions.R')
    source('VarianceFunctions.R')
    source('NonStGM_adj.R')
    source('parallel_sim_dependencies.R')
    library(dplyr)
    library(tibble)
  })

  clusterExport(cl, c("sim.tvVAR_p3"), envir = environment())
  clusterSetRNGStream(cl, seed)

  All_Test_tibble = foreach(iter = 1:R, .combine = rbind, .packages = c("dplyr", "tibble")) %dopar% {
    x = sim.tvVAR_p3(burnin, m, TV_size)
    tib = NonStGM_adj(x, Kernel = Kernel, nu = nu, L = L, alpha = alpha)
    tib = tib |> mutate(i = iter)
    tib
  }

  stopCluster(cl)

  results = graph_recovery_p3(All_Test_tibble, alpha = alpha)
  return(list(results, All_Test_tibble))
}

# ============================================================
# RUN THE SIMULATION (PARALLEL)
# ============================================================
results_p3 = run_NonStGM_simulation_p3_parallel(
  R = 1000,
  alpha = 0.05,
  burnin = 200,
  m = 2^12,
  TV_size = 0.6,
  nu = 2,
  L = 1,
  Kernel = 'Kernel_Triangular',
  seed = 1,
  n_cores = 8
)

print(results_p3[[1]])

# ============================================================
# PARALLEL SIMULATION RUNNER FOR p=5
# ============================================================
run_NonStGM_simulation_p5_parallel = function(R = 100, alpha = 0.05, burnin = 20, m = 4096,
                                              TV_size = 0.6, nu = 2, L = 1,
                                              Kernel = 'Kernel_Triangular', seed = 0,
                                              n_cores = NULL) {

  if (is.null(n_cores)) {
    n_cores = detectCores() - 1
  }

  cl = makeCluster(n_cores)
  registerDoParallel(cl)

  # Necessary dependencies
  clusterEvalQ(cl, {
    source('extractVariables.R')
    source('exampleVAR.R')
    source('KernelWeights.R')
    source('Combined_MK_Estimation.R')
    source('DFTransform.R')
    source('R_Hat_Creation.R')
    source('BetaFunctions.R')
    source('VarianceFunctions.R')
    source('NonStGM_adj.R')
    source('parallel_sim_dependencies.R')
    library(dplyr)
    library(tibble)
  })

  # Export the simulation function
  clusterExport(cl, c("sim.tvVAR_p5"), envir = environment())

  clusterSetRNGStream(cl, seed)

  All_Test_tibble = foreach(iter = 1:R, .combine = rbind, .packages = c("dplyr", "tibble")) %dopar% {
    x = sim.tvVAR_p5(burnin, m, TV_size)
    tib = NonStGM_adj(x, Kernel = Kernel, nu = nu, L = L, alpha = alpha)
    tib = tib |> mutate(i = iter)
    tib
  }

  stopCluster(cl)

  results = graph_recovery_p5(All_Test_tibble, alpha = alpha)
  return(results)
}

# ============================================================
# PARALLEL SIMULATION RUNNER FOR p=7
# ============================================================
run_NonStGM_simulation_p7_parallel = function(R = 100, alpha = 0.05, burnin = 20, m = 4096,
                                              TV_size = 0.6, nu = 2, L = 1,
                                              Kernel = 'Kernel_Triangular', seed = 0,
                                              n_cores = NULL) {

  if (is.null(n_cores)) {
    n_cores = detectCores() - 1
  }

  cl = makeCluster(n_cores)
  registerDoParallel(cl)

  clusterEvalQ(cl, {
    source('extractVariables.R')
    source('exampleVAR.R')
    source('KernelWeights.R')
    source('Combined_MK_Estimation.R')
    source('DFTransform.R')
    source('R_Hat_Creation.R')
    source('BetaFunctions.R')
    source('VarianceFunctions.R')
    source('NonStGM_adj.R')
    source('parallel_sim_dependencies.R')
    library(dplyr)
    library(tibble)
  })

  clusterExport(cl, c("sim.tvVAR_p7"), envir = environment())
  clusterSetRNGStream(cl, seed)

  All_Test_tibble = foreach(iter = 1:R, .combine = rbind, .packages = c("dplyr", "tibble")) %dopar% {
    x = sim.tvVAR_p7(burnin, m, TV_size)
    tib = NonStGM_adj(x, Kernel = Kernel, nu = nu, L = L, alpha = alpha)
    tib = tib |> mutate(i = iter)
    tib
  }

  stopCluster(cl)

  results = graph_recovery_p7(All_Test_tibble, alpha = alpha)
  return(list(results, All_Test_tibble))
}

# ============================================================
# RUN PARALLEL SIMULATIONS
# ============================================================
results_p5 = run_NonStGM_simulation_p5_parallel(
  R = 100,
  alpha = 0.05,
  burnin = 200,
  m = 2^12,
  TV_size = 0.6,
  nu = 2,
  L = 1,
  Kernel = 'Kernel_Triangular',
  seed = 1,
  n_cores = 8
)

print(results_p5)

results_p7 = run_NonStGM_simulation_p7_parallel(
  R = 100,
  alpha = 0.05,
  burnin = 200,
  m = 2^12,
  TV_size = .6,
  nu = 2,
  L = 1,
  Kernel = 'Kernel_Triangular',
  seed = 1,
  n_cores = 8
)
print(results_p7[[1]])

# filtered_results = results_p7[[2]] |>
#   filter(a == 1 & c == 1 & r == 1)
#
# p_value_distribution = filtered_results |>
#   group_by(a,c,i) |>
#   summarise(
#     p_min = min(p.adjust(c(Re, Im), method = "BY")),
#     .groups = "drop"
#   )
#
# # View distribution
# p_list = p_value_distribution
#
# # Histogram
# hist(p_value_distribution$p_min, breaks = 20,
#      main = "P-value Distribution (a=1, c=1)",
#      xlab = "BY-adjusted p-value")




run_NonStGM_simulation2_p7_parallel = function(R = 100, alpha = 0.05, burnin = 20, m = 4096,
                                               TV_size = 0.6, nu = 2, L = 1,
                                               Kernel = 'Kernel_Triangular', seed = 0,
                                               n_cores = NULL) {
  if (is.null(n_cores)) {
    n_cores = detectCores() - 1
  }
  cl = makeCluster(n_cores)
  registerDoParallel(cl)
  clusterEvalQ(cl, {
    source('extractVariables.R')
    source('exampleVAR.R')
    source('KernelWeights.R')
    source('Combined_MK_Estimation.R')
    source('DFTransform.R')
    source('R_Hat_Creation.R')
    source('BetaFunctions.R')
    source('VarianceFunctions.R')
    source('NonStGM_adj.R')
    source('parallel_sim_dependencies.R')
    library(dplyr)
    library(tibble)
  })
  clusterExport(cl, c("sim.tvVAR2_p7"), envir = environment())
  clusterSetRNGStream(cl, seed)
  All_Test_tibble = foreach(iter = 1:R, .combine = rbind, .packages = c("dplyr", "tibble")) %dopar% {
    x = sim.tvVAR2_p7(burnin, m, TV_size)
    tib = NonStGM_adj(x, Kernel = Kernel, nu = nu, L = L, alpha = alpha)
    tib = tib |> mutate(i = iter)
    tib
  }
  stopCluster(cl)
  results = graph_recovery2_p7(All_Test_tibble, alpha = alpha)
  return(list(results, All_Test_tibble))
}

results2_p7 = run_NonStGM_simulation2_p7_parallel(
  R = 1000,
  alpha = 0.05,
  burnin = 200,
  m = 2^12,
  TV_size = .6,
  nu = 2,
  L = 1,
  Kernel = 'Kernel_Triangular',
  seed = 1,
  n_cores = 8
)

print(results2_p7[[1]])

# run_NonStGM_simulation_p15_parallel = function(R = 100, alpha = 0.05, burnin = 20, m = 4096,
#                                                TV_size = 0.6, nu = 2, L = 1,
#                                                Kernel = 'Kernel_Triangular', seed = 0,
#                                                n_cores = NULL) {
#   if (is.null(n_cores)) {
#     n_cores = detectCores() - 1
#   }
#   cl = makeCluster(n_cores)
#   registerDoParallel(cl)
#   clusterEvalQ(cl, {
#     source('extractVariables.R')
#     source('exampleVAR.R')
#     source('KernelWeights.R')
#     source('Combined_MK_Estimation.R')
#     source('DFTransform.R')
#     source('R_Hat_Creation.R')
#     source('BetaFunctions.R')
#     source('VarianceFunctions.R')
#     source('NonStGM_adj.R')
#     source('parallel_sim_dependencies.R')
#     library(dplyr)
#     library(tibble)
#   })
#   clusterExport(cl, c("sim.tvVAR_p15"), envir = environment())
#   clusterSetRNGStream(cl, seed)
#   All_Test_tibble = foreach(iter = 1:R, .combine = rbind, .packages = c("dplyr", "tibble")) %dopar% {
#     x = sim.tvVAR_p15(burnin, m, TV_size)
#     tib = NonStGM_adj(x, Kernel = Kernel, nu = nu, L = L, alpha = alpha)
#     tib = tib |> mutate(i = iter)
#     tib
#   }
#   stopCluster(cl)
#   results = graph_recovery_p15(All_Test_tibble, alpha = alpha)
#   return(list(results, All_Test_tibble))
# }
#
# results_p15 = run_NonStGM_simulation_p15_parallel(
#   R = 10,
#   alpha = 0.05,
#   burnin = 200,
#   m = 2^12,
#   TV_size = .6,
#   nu = 2,
#   L = 1,
#   Kernel = 'Kernel_Triangular',
#   seed = 1,
#   n_cores = 8
# )
# print(results_p15[[1]])
