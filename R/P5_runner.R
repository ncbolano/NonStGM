library(foreach)
library(doParallel)
library(parallel)

source('full_function_list.R')

p = 5


run_NonStGM_simulation_p5_parallel_compare = function(R = 100, alpha = 0.05, burnin = 20, m = 4096,
                                                      TV_size = 0.6, nu = 2, L = 1,
                                                      Kernel = 'Kernel_Triangular', seed = 0,
                                                      n_cores = NULL) {
  if (is.null(n_cores)) n_cores = detectCores() - 1
  cl = makeCluster(n_cores)
  registerDoParallel(cl)

  clusterEvalQ(cl, {
    source('full_function_list.R')
    library(dplyr)
    library(tibble)
  })

  clusterExport(cl, c(
    "beta", "beta_beta_r", "beta_cv2",
    "extractVariables", "extractK", "extractM", "local_M_selection",
    "return_max", "r_to_loc",
    "Kernel_Triangular", "Kernel_Quadratic", "Kernel_Uniform",
    "Hat.r", "Hat.R", "Hat.R.reduced", "J.to.J.k.nu", "J.to.J.k.nu.a",
    "variance.estimator.v2", "variance.estimator.Re.v2", "variance.estimator.Im.v2",
    "whittle_loss", "smoothed_spectral_density", "smoothed_spectral_density_LOO",
    "sim.tvVAR_p5", "NonStGM_adj", "extract_skeleton", "NonStGM_node_refit",
    "graph_recovery_p5", "graph_recovery_p5_twostep"
  ), envir = environment())

  clusterSetRNGStream(cl, seed)

  combine_lists = function(x, y) list(
    full  = rbind(x$full,  y$full),
    edges = rbind(x$edges, y$edges),
    nodes = rbind(x$nodes, y$nodes)
  )

  combined = foreach(iter = 1:R, .combine = combine_lists,
                     .init = list(full = NULL, edges = NULL, nodes = NULL),
                     .packages = c("dplyr", "tibble")) %dopar% {
                       x = sim.tvVAR_p5(burnin, m, TV_size)
                       tib_full = NonStGM_adj(x, Kernel = Kernel, nu = nu, L = L, alpha = alpha,
                                              skeleton_only = FALSE) |> mutate(i = iter)
                       tib_edge = tib_full |> filter(r == 0)
                       nb = extract_skeleton(tib_edge, p = 5, alpha = alpha)
                       tib_node = NonStGM_node_refit(x, neighbors = nb, Kernel = Kernel,
                                                     nu = nu, L = L, alpha = alpha) |> mutate(i = iter)
                       list(full = tib_full, edges = tib_edge, nodes = tib_node)
                     }
  stopCluster(cl)

  list(
    onestep = graph_recovery_p5(combined$full, alpha = alpha),
    twostep = graph_recovery_p5_twostep(combined$edges, combined$nodes, alpha = alpha),
    full_tibble = combined$full,
    edge_tibble = combined$edges,
    node_tibble = combined$nodes
  )
}

results_p5_compare = run_NonStGM_simulation_p5_parallel_compare(
  R = 500, alpha = 0.01, burnin = 200, m = 2^12,
  TV_size = 0.45, nu = 2, L = 1,
  Kernel = 'Kernel_Triangular', seed = 1, n_cores = 20
)

print(bind_rows(results_p5_compare$onestep, results_p5_compare$twostep))




















################ Test 2 (all time varying nodes)

run_NonStGM_simulation_p5_alltv_parallel_compare = function(R = 100, alpha = 0.05, burnin = 20, m = 4096,
                                                      TV_size = 0.6, nu = 2, L = 1,
                                                      Kernel = 'Kernel_Triangular', seed = 0,
                                                      n_cores = NULL) {
  if (is.null(n_cores)) n_cores = detectCores() - 1
  cl = makeCluster(n_cores)
  registerDoParallel(cl)

  clusterEvalQ(cl, {
    source('full_function_list.R')
    library(dplyr)
    library(tibble)
  })

  clusterExport(cl, c(
    "beta", "beta_beta_r", "beta_cv2",
    "extractVariables", "extractK", "extractM", "local_M_selection",
    "return_max", "r_to_loc",
    "Kernel_Triangular", "Kernel_Quadratic", "Kernel_Uniform",
    "Hat.r", "Hat.R", "Hat.R.reduced", "J.to.J.k.nu", "J.to.J.k.nu.a",
    "variance.estimator.v2", "variance.estimator.Re.v2", "variance.estimator.Im.v2",
    "whittle_loss", "smoothed_spectral_density", "smoothed_spectral_density_LOO",
    "sim.tvVAR_p5_alltv", "NonStGM_adj", "extract_skeleton", "NonStGM_node_refit",
    "graph_recovery_p5_alltv", "graph_recovery_p5_alltv_twostep"
  ), envir = environment())

  clusterSetRNGStream(cl, seed)

  combine_lists = function(x, y) list(
    full  = rbind(x$full,  y$full),
    edges = rbind(x$edges, y$edges),
    nodes = rbind(x$nodes, y$nodes)
  )

  combined = foreach(iter = 1:R, .combine = combine_lists,
                     .init = list(full = NULL, edges = NULL, nodes = NULL),
                     .packages = c("dplyr", "tibble")) %do% {
                       cat("simulation\n")
                       x = sim.tvVAR_p5_alltv(burnin, m, TV_size)
                       cat("NonStGM_adj\n")
                       tib_full = NonStGM_adj(x, Kernel = Kernel, nu = nu, L = L, alpha = alpha,
                                              skeleton_only = FALSE) |> mutate(i = iter)
                       tib_edge = tib_full |> filter(r == 0)
                       cat("extract skeleton\n")
                       nb = extract_skeleton(tib_edge, p = 5, alpha = alpha)
                       cat("node refit\n")
                       tib_node = NonStGM_node_refit(x, neighbors = nb, Kernel = Kernel,
                                                     nu = nu, L = L, alpha = alpha) |> mutate(i = iter)
                       list(full = tib_full, edges = tib_edge, nodes = tib_node)
                     }
  stopCluster(cl)

  list(
    onestep = graph_recovery_p5_alltv(combined$full, alpha = alpha),
    twostep = graph_recovery_p5_alltv_twostep(combined$edges, combined$nodes, alpha = alpha),
    full_tibble = combined$full,
    edge_tibble = combined$edges,
    node_tibble = combined$nodes
  )
}

results_p5_alltv_compare = run_NonStGM_simulation_p5_alltv_parallel_compare(
  R = 500, alpha = 0.01, burnin = 200, m = 2^12,
  TV_size = 0.3, nu = 2, L = 1,
  Kernel = 'Kernel_Triangular', seed = 1, n_cores = 20
)


print(bind_rows(results_p5_alltv_compare$onestep, results_p5_alltv_compare$twostep))
