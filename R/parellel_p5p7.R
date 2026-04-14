library(future)
library(furrr)

library(foreach)
library(doParallel)
library(parallel)

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
  n_cores = 6
)
print(results_p5)

results_p7 = run_NonStGM_simulation_p7_parallel(
  R = 100,
  alpha = 0.05,
  burnin = 20,
  m = 2^12,
  TV_size = 0.6,
  nu = 2,
  L = 1,
  Kernel = 'Kernel_Triangular',
  seed = 0,
  n_cores = 6
)
print(results_p7[[1]])
