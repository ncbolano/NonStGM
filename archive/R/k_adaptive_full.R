
source('extractVariables.R')
source('exampleVAR.R')
source('KernelWeights.R')
source('Combined_MK_Estimation.R')
source('DFTransform.R')
source('R_Hat_Creation.R')
source('BetaFunctions.R')
source('VarianceFunctions.R')

library(foreach)
library(doParallel)
library(parallel)
library(tidyverse)

p = 7


sim.tvVAR_p7 = function(burnin, m, TV_size) {
  A = matrix(c(
    0.5,  0,    0,    0,    0,    0,    0,
    0.3,  0.5,  0,    0,    0,    0,    0,
    0,    0.3,  0.5,  0,    0,    0,    0,
    0,    0,    0.3,  0.5,  0,    0,    0,
    0,    0,    0,    0.3,  0.5,  0,    0,
    0,    0,    0,    0,    0.3,  0.5,  0,
    0.3,  0,    0,    0,    0,    0.3,  0.5
  ), ncol = 7, byrow = TRUE)
  n = m + burnin
  p = 7
  x = matrix(rnorm(n * p), ncol = p)
  x1 = x
  st = TV_size * (1 + exp(0.005 * (c(1:n) - (n / 2))))^(-1)
  for (tt in 2:n) {
    A.t = A
    A.t[1, 1] = st[tt]
    temp = A.t %*% matrix(x1[tt - 1, ], ncol = 1) + matrix(x[tt, ], ncol = 1)
    x1[tt, ] = c(temp)
  }
  x2 = x1[-c(1:burnin), ]
  return(x2)
}


NonStGM_edges_only = function(x, Kernel = 'Kernel_Triangular', nu = 2, L = 1, alpha = 0.05) {
  Kernel_function = get(Kernel)
  variable_list = extractVariables(x)
  J = variable_list[[1]]
  p = as.numeric(variable_list[2])
  coefnum = (2 * p * nu) + p - 1
  k = extractK(J, coefnum)
  n = nrow(x)
  M_grid = unique(round(seq(max(coefnum, n^(1/5)), max(coefnum, n^(1/2)), length.out = 10)))
  M_list = local_M_selection(J, k, M_grid)
  n_k = length(k)
  betaCoefAll = array(0, c(n_k, coefnum, p))
  varbeta = array(0, c(n_k, coefnum, 2, p))
  varbetaAll = array(0, c(n_k, 2 * coefnum, 2 * coefnum, p))
  for (j in seq_along(k)) {
    for (i_p in 1:p) {
      betaCoefAll[j, , i_p] = beta(J, k[j], nu, Kernel_function, M_list[j], a = i_p, delta = 0)
      varbetaAll[j, , , i_p] = variance.estimator.v2(J, k[j], nu, Kernel_function, M_list[j], a = i_p, L, delta = 0)
      varbeta[j, , 1, i_p] = diag(varbetaAll[j, 1:coefnum, 1:coefnum, i_p])
      varbeta[j, , 2, i_p] = diag(varbetaAll[j, 1:coefnum + coefnum, 1:coefnum + coefnum, i_p])
    }
  }
  Test_tibble = NULL
  for (a in 1:p) {
    for (c in 1:p) {
      if (a < c) {
        r = 0
        for (j in seq_along(k)) {
          loc = r_to_loc(r, c, a, nu = nu, p = p)
          tmp1 = pnorm(abs(Re(betaCoefAll[j, loc, a])) / sqrt(abs(varbeta[j, loc, 1, a])), lower.tail = F) * 2
          tmp2 = pnorm(abs(Im(betaCoefAll[j, loc, a])) / sqrt(abs(varbeta[j, loc, 2, a])), lower.tail = F) * 2
          Test_tibble = tibble(Re = tmp1, Im = tmp2, a = a, c = c, r = r, k = k[j], i = 1) |>
            rbind(Test_tibble)
        }
      }
    }
  }
  return(list(Test_tibble = Test_tibble, k = k, M_list = M_list))
}

edge_recovery_p7 = function(Test_tibble, alpha = 0.05) {
  Graph_tibble = Test_tibble |>
    group_by(a, c, i) |>
    summarise(
      p_min = min(p.adjust(c(Re, Im), method = "BY")),
      .groups = "drop"
    ) |>
    mutate(edge = as.integer(p_min < alpha))
  Edge_list = Graph_tibble |>
    filter(edge == 1) |>
    group_by(i) |>
    summarise(
      neighbors = list(tibble(a = a, c = c)),
      .groups = "drop"
    )
  edge_results = Graph_tibble |>
    group_by(i) |>
    summarise(
      detect_12 = any(a == 1 & c == 2 & p_min < alpha),
      detect_23 = any(a == 2 & c == 3 & p_min < alpha),
      detect_34 = any(a == 3 & c == 4 & p_min < alpha),
      detect_45 = any(a == 4 & c == 5 & p_min < alpha),
      detect_56 = any(a == 5 & c == 6 & p_min < alpha),
      detect_67 = any(a == 6 & c == 7 & p_min < alpha),
      detect_17 = any(a == 1 & c == 7 & p_min < alpha),
      n_offdiag_edges = sum(p_min < alpha),
      .groups = "drop"
    )
  edge_accuracy = tibble(
    Detect_12 = mean(edge_results$detect_12),
    Detect_23 = mean(edge_results$detect_23),
    Detect_34 = mean(edge_results$detect_34),
    Detect_45 = mean(edge_results$detect_45),
    Detect_56 = mean(edge_results$detect_56),
    Detect_67 = mean(edge_results$detect_67),
    Detect_17 = mean(edge_results$detect_17),
    Mean_N_Edges = mean(edge_results$n_offdiag_edges),
    Mean_False_Edges = mean(pmax(0, edge_results$n_offdiag_edges - 7))
  )
  return(list(
    edge_accuracy = edge_accuracy,
    Graph_tibble = Graph_tibble,
    Edge_list = Edge_list
  ))
}


test_node_tv = function(x, node, neighbors, k_full, M_list_full,
                        Kernel = 'Kernel_Triangular', nu = 2, L = 1, alpha = 0.05) {
  subset_cols = sort(unique(c(node, neighbors)))
  x_sub = x[, subset_cols]
  node_sub = which(subset_cols == node)
  Kernel_function = get(Kernel)
  variable_list = extractVariables(x_sub)
  J = variable_list[[1]]
  p_sub = as.numeric(variable_list[2])
  coefnum = (2 * p_sub * nu) + p_sub - 1
  k = k_full
  M_list = M_list_full
  n_k = length(k)
  betaCoefAll = array(0, c(n_k, coefnum, p_sub))
  varbeta = array(0, c(n_k, coefnum, 2, p_sub))
  varbetaAll = array(0, c(n_k, 2 * coefnum, 2 * coefnum, p_sub))
  for (j in seq_along(k)) {
    betaCoefAll[j, , node_sub] = beta(J, k[j], nu, Kernel_function, M_list[j], a = node_sub, delta = 0)
    varbetaAll[j, , , node_sub] = variance.estimator.v2(J, k[j], nu, Kernel_function, M_list[j], a = node_sub, L, delta = 0)
    varbeta[j, , 1, node_sub] = diag(varbetaAll[j, 1:coefnum, 1:coefnum, node_sub])
    varbeta[j, , 2, node_sub] = diag(varbetaAll[j, 1:coefnum + coefnum, 1:coefnum + coefnum, node_sub])
  }
  TV_tibble = NULL
  a = node_sub
  c_idx = node_sub
  for (r in 1:nu) {
    for (j in seq_along(k)) {
      loc = r_to_loc(r, c_idx, a, nu = nu, p = p_sub)
      tmp1 = pnorm(abs(Re(betaCoefAll[j, loc, a])) / sqrt(abs(varbeta[j, loc, 1, a])), lower.tail = F) * 2
      tmp2 = pnorm(abs(Im(betaCoefAll[j, loc, a])) / sqrt(abs(varbeta[j, loc, 2, a])), lower.tail = F) * 2
      TV_tibble = tibble(Re = tmp1, Im = tmp2, node = node, r = r, k = k[j], p_sub = p_sub) |>
        rbind(TV_tibble)
    }
  }
  return(TV_tibble)
}

twostep_recovery_p7 = function(step1_results, TV_results, alpha = 0.05) {
  edge_graph = step1_results$Graph_tibble
  tv_results = TV_results
  edge_summary = edge_graph |>
    group_by(i) |>
    summarise(
      detect_12 = any(a == 1 & c == 2 & p_min < alpha),
      detect_23 = any(a == 2 & c == 3 & p_min < alpha),
      detect_34 = any(a == 3 & c == 4 & p_min < alpha),
      detect_45 = any(a == 4 & c == 5 & p_min < alpha),
      detect_56 = any(a == 5 & c == 6 & p_min < alpha),
      detect_67 = any(a == 6 & c == 7 & p_min < alpha),
      detect_17 = any(a == 1 & c == 7 & p_min < alpha),
      n_offdiag_edges = sum(p_min < alpha),
      .groups = "drop"
    )
  tv_summary = tv_results |>
    group_by(i) |>
    summarise(
      detect_11_tv = any(node == 1 & tv_detected),
      false_tv_2 = any(node == 2 & tv_detected),
      false_tv_3 = any(node == 3 & tv_detected),
      false_tv_4 = any(node == 4 & tv_detected),
      false_tv_5 = any(node == 5 & tv_detected),
      false_tv_6 = any(node == 6 & tv_detected),
      false_tv_7 = any(node == 7 & tv_detected),
      n_tv_detected = sum(tv_detected),
      .groups = "drop"
    )
  combined = edge_summary |>
    left_join(tv_summary, by = "i")
  accuracy_summary = combined |>
    mutate(
      all_offdiag_detected = detect_12 & detect_23 & detect_34 &
        detect_45 & detect_56 & detect_67 & detect_17,
      all_true_detected = all_offdiag_detected & detect_11_tv,
      no_false_edges = (n_offdiag_edges == 7),
      no_false_tv = (n_tv_detected <= 1),
      perfect = all_offdiag_detected & no_false_edges & detect_11_tv & no_false_tv,
      sensitivity = (detect_11_tv + detect_12 + detect_23 + detect_34 +
                       detect_45 + detect_56 + detect_67 + detect_17) / 8,
      false_edges = pmax(0, n_offdiag_edges - 7),
      false_tv = pmax(0, n_tv_detected - 1)
    )
  edge_rates = c(
    mean(accuracy_summary$detect_12),
    mean(accuracy_summary$detect_23),
    mean(accuracy_summary$detect_34),
    mean(accuracy_summary$detect_45),
    mean(accuracy_summary$detect_56),
    mean(accuracy_summary$detect_67),
    mean(accuracy_summary$detect_17)
  )
  overall_accuracy = tibble(
    Method = "Two-Step",
    Perfect_Recovery = mean(accuracy_summary$perfect),
    All_True_Detected = mean(accuracy_summary$all_true_detected),
    Smallest_Edge = min(edge_rates),
    Detect_11_TV = mean(accuracy_summary$detect_11_tv),
    Mean_Sensitivity = mean(accuracy_summary$sensitivity),
    Detect_12 = mean(accuracy_summary$detect_12),
    Detect_23 = mean(accuracy_summary$detect_23),
    Detect_34 = mean(accuracy_summary$detect_34),
    Detect_45 = mean(accuracy_summary$detect_45),
    Detect_56 = mean(accuracy_summary$detect_56),
    Detect_67 = mean(accuracy_summary$detect_67),
    Detect_17 = mean(accuracy_summary$detect_17),
    False_TV_rate = mean(accuracy_summary$false_tv > 0),
    Mean_False_TV = mean(accuracy_summary$false_tv),
    Mean_N_Edges = mean(accuracy_summary$n_offdiag_edges),
    Mean_False_Edges = mean(accuracy_summary$false_edges),
    Mean_Total_FP = mean(accuracy_summary$false_edges + accuracy_summary$false_tv)
  )
  return(overall_accuracy)
}

run_twostep_p7 = function(R = 100, alpha = 0.05, burnin = 200, m = 4096,
                          TV_size = 0.6, nu = 2, L = 1,
                          Kernel = 'Kernel_Triangular', seed = 1,
                          n_cores = 4) {

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
    library(dplyr)
    library(tibble)
  })
  clusterExport(cl, c("sim.tvVAR_p7", "NonStGM_edges_only", "r_to_loc"), envir = environment())
  clusterSetRNGStream(cl, seed)

  all_step1 = foreach(iter = 1:R, .packages = c("dplyr", "tibble")) %dopar% {
    x = sim.tvVAR_p7(burnin, m, TV_size)
    res = NonStGM_edges_only(x, Kernel = Kernel, nu = nu, L = L, alpha = alpha)
    tib = res$Test_tibble |> mutate(i = iter)
    list(tibble = tib, k = res$k, M_list = res$M_list, x = x)
  }
  stopCluster(cl)

  All_Test_tibble = bind_rows(lapply(all_step1, function(r) r$tibble))
  k_list = lapply(all_step1, function(r) r$k)
  M_list_list = lapply(all_step1, function(r) r$M_list)
  x_list = lapply(all_step1, function(r) r$x)

  step1 = edge_recovery_p7(All_Test_tibble, alpha = alpha)
  print(step1$edge_accuracy)

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
    library(dplyr)
    library(tibble)
  })
  clusterExport(cl, c("test_node_tv", "r_to_loc"), envir = environment())
  clusterSetRNGStream(cl, 42)

  Edge_list = step1$Edge_list

  tv_tibble = foreach(iter = 1:R, .combine = rbind, .packages = c("dplyr", "tibble")) %dopar% {
    x = x_list[[iter]]
    k_full = k_list[[iter]]
    M_list_full = M_list_list[[iter]]
    edges_i = Edge_list |> filter(i == iter) |> pull(neighbors)
    if (length(edges_i) == 0) return(NULL)
    edges_df = edges_i[[1]]
    iter_results = NULL
    for (node in 1:7) {
      neighbors = unique(c(
        edges_df$c[edges_df$a == node],
        edges_df$a[edges_df$c == node]
      ))
      if (length(neighbors) == 0) next
      tv_tib = test_node_tv(x, node, neighbors, k_full, M_list_full,
                            Kernel = Kernel, nu = nu, L = L, alpha = alpha)
      tv_tib = tv_tib |> mutate(i = iter)
      iter_results = rbind(iter_results, tv_tib)
    }
    iter_results
  }
  stopCluster(cl)

  TV_results = tv_tibble |>
    group_by(node, i) |>
    summarise(
      p_min = min(p.adjust(c(Re, Im), method = "BY")),
      p_sub_used = first(p_sub),
      .groups = "drop"
    ) |>
    mutate(tv_detected = p_min < alpha)

  TV_summary = TV_results |>
    group_by(node) |>
    summarise(
      detection_rate = mean(tv_detected),
      mean_p_sub = mean(p_sub_used),
      .groups = "drop"
    )

  cat("\nStep 2 TV Summary:\n")
  print(TV_summary)

  final = twostep_recovery_p7(step1, TV_results, alpha = alpha)
  cat("\nFinal Results:\n")
  print(final)

  return(list(
    step1 = step1,
    tv_tibble = tv_tibble,
    TV_results = TV_results,
    TV_summary = TV_summary,
    final = final
  ))
}


results = run_twostep_p7(
  R = 100,
  alpha = 0.05,
  burnin = 200,
  m = 2^12,
  TV_size = 0.6,
  nu = 2,
  L = 1,
  Kernel = 'Kernel_Triangular',
  seed = 1,
  n_cores = 16
)
