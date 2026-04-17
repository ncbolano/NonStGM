twostep_recovery_p7 = function(step1_results, TV_results, alpha = 0.05) {

  edge_graph = step1_results$Graph_tibble  # i, a, c, p_min, edge
  tv_results = TV_results                   # node, i, p_min, tv_detected

  # Per-replication edge summary
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

  # Per-replication TV summary (only node 1 is truly TV)
  tv_summary = tv_results |>
    group_by(i) |>
    summarise(
      detect_11_tv = any(node == 1 & tv_detected),
      # False TV detections (nodes 2-7 should NOT be TV)
      false_tv_2 = any(node == 2 & tv_detected),
      false_tv_3 = any(node == 3 & tv_detected),
      false_tv_4 = any(node == 4 & tv_detected),
      false_tv_5 = any(node == 5 & tv_detected),
      false_tv_6 = any(node == 6 & tv_detected),
      false_tv_7 = any(node == 7 & tv_detected),
      n_tv_detected = sum(tv_detected),
      .groups = "drop"
    )

  # Join edge + TV results
  combined = edge_summary |>
    left_join(tv_summary, by = "i")

  # Accuracy metrics (SAME as your original table)
  accuracy_summary = combined |>
    mutate(
      all_offdiag_detected = detect_12 & detect_23 & detect_34 &
        detect_45 & detect_56 & detect_67 & detect_17,
      all_true_detected = all_offdiag_detected & detect_11_tv,
      no_false_edges = (n_offdiag_edges == 7),
      no_false_tv = (n_tv_detected <= 1),  # Only node 1 should be TV
      perfect = all_offdiag_detected & no_false_edges & detect_11_tv & no_false_tv,
      sensitivity = (detect_11_tv + detect_12 + detect_23 + detect_34 +
                       detect_45 + detect_56 + detect_67 + detect_17) / 8,
      false_edges = pmax(0, n_offdiag_edges - 7),
      false_tv = pmax(0, n_tv_detected - 1)
    )

  # Smallest edge detection
  edge_rates = c(
    mean(accuracy_summary$detect_12),
    mean(accuracy_summary$detect_23),
    mean(accuracy_summary$detect_34),
    mean(accuracy_summary$detect_45),
    mean(accuracy_summary$detect_56),
    mean(accuracy_summary$detect_67),
    mean(accuracy_summary$detect_17)
  )

  # Overall summary (MATCHES YOUR TABLE FORMAT)
  overall_accuracy = tibble(
    Method = "Two-Step",
    Perfect_Recovery = mean(accuracy_summary$perfect),
    All_True_Detected = mean(accuracy_summary$all_true_detected),
    Smallest_Edge = min(edge_rates),
    Detect_11_TV = mean(accuracy_summary$detect_11_tv),
    Mean_Sensitivity = mean(accuracy_summary$sensitivity),
    # Edge details
    Detect_12 = mean(accuracy_summary$detect_12),
    Detect_23 = mean(accuracy_summary$detect_23),
    Detect_34 = mean(accuracy_summary$detect_34),
    Detect_45 = mean(accuracy_summary$detect_45),
    Detect_56 = mean(accuracy_summary$detect_56),
    Detect_67 = mean(accuracy_summary$detect_67),
    Detect_17 = mean(accuracy_summary$detect_17),
    # TV details
    False_TV_rate = mean(accuracy_summary$false_tv > 0),
    Mean_False_TV = mean(accuracy_summary$false_tv),
    # Edge FP
    Mean_N_Edges = mean(accuracy_summary$n_offdiag_edges),
    Mean_False_Edges = mean(accuracy_summary$false_edges),
    # Combined FP
    Mean_Total_FP = mean(accuracy_summary$false_edges + accuracy_summary$false_tv)
  )

  return(overall_accuracy)
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
          Test_tibble = tibble(
            Re = tmp1,
            Im = tmp2,
            a = a,
            c = c,
            r = r,
            k = k[j],
            i = 1
          ) |> rbind(Test_tibble)
        }
      }
    }
  }

  # RETURN k and M_list alongside tibble
  return(list(Test_tibble = Test_tibble, k = k, M_list = M_list))
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

  # USE STEP 1 k AND M VALUES DIRECTLY
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
      TV_tibble = tibble(
        Re = tmp1,
        Im = tmp2,
        node = node,
        r = r,
        k = k[j],
        p_sub = p_sub
      ) |> rbind(TV_tibble)
    }
  }

  return(TV_tibble)
}

run_edge_detection_p7_parallel = function(R = 100, alpha = 0.05, burnin = 20, m = 4096,
                                          TV_size = 0.6, nu = 2, L = 1,
                                          Kernel = 'Kernel_Triangular', seed = 1,
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
    library(dplyr)
    library(tibble)
  })

  clusterExport(cl, c("sim.tvVAR_p7", "NonStGM_edges_only", "r_to_loc"), envir = environment())
  clusterSetRNGStream(cl, seed)

  # Return list instead of rbind — need to keep k and M per iteration
  all_results = foreach(iter = 1:R, .packages = c("dplyr", "tibble")) %dopar% {
    x = sim.tvVAR_p7(burnin, m, TV_size)
    res = NonStGM_edges_only(x, Kernel = Kernel, nu = nu, L = L, alpha = alpha)
    tib = res$Test_tibble |> mutate(i = iter)
    list(tibble = tib, k = res$k, M_list = res$M_list, x = x)
  }

  stopCluster(cl)

  # Extract components
  All_Test_tibble = bind_rows(lapply(all_results, function(r) r$tibble))
  k_list = lapply(all_results, function(r) r$k)
  M_list_list = lapply(all_results, function(r) r$M_list)
  x_list = lapply(all_results, function(r) r$x)

  step1 = edge_recovery_p7(All_Test_tibble, alpha = alpha)

  return(list(
    edge_accuracy = step1$edge_accuracy,
    Graph_tibble = step1$Graph_tibble,
    Edge_list = step1$Edge_list,
    raw_tibble = All_Test_tibble,
    k_list = k_list,
    M_list_list = M_list_list,
    x_list = x_list
  ))
}

run_node_tv_p7_parallel = function(x_list, Edge_list, k_list, M_list_list,
                                   nodes_to_test = 1:7,
                                   Kernel = 'Kernel_Triangular', nu = 2, L = 1,
                                   alpha = 0.05, n_cores = NULL) {
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
    library(dplyr)
    library(tibble)
  })

  clusterExport(cl, c("test_node_tv", "r_to_loc"), envir = environment())
  clusterSetRNGStream(cl, 42)

  R = length(x_list)

  All_TV_tibble = foreach(iter = 1:R, .combine = rbind, .packages = c("dplyr", "tibble")) %dopar% {

    x = x_list[[iter]]
    k_full = k_list[[iter]]
    M_list_full = M_list_list[[iter]]

    edges_i = Edge_list |> filter(i == iter) |> pull(neighbors)
    if (length(edges_i) == 0) return(NULL)
    edges_df = edges_i[[1]]

    iter_results = NULL
    for (node in nodes_to_test) {
      neighbors = unique(c(
        edges_df$c[edges_df$a == node],
        edges_df$a[edges_df$c == node]
      ))

      if (length(neighbors) == 0) next

      # PASS k and M from Step 1
      tv_tib = test_node_tv(x, node, neighbors, k_full, M_list_full,
                            Kernel = Kernel, nu = nu, L = L, alpha = alpha)
      tv_tib = tv_tib |> mutate(i = iter)
      iter_results = rbind(iter_results, tv_tib)
    }
    iter_results
  }

  stopCluster(cl)
  return(All_TV_tibble)
}

run_twostep_p7 = function(R = 100, alpha = 0.05, burnin = 200, m = 4096,
                          TV_size = 0.6, nu = 2, L = 1,
                          Kernel = 'Kernel_Triangular', seed = 1,
                          n_cores = 20) {

  step1 = run_edge_detection_p7_parallel(
    R = R, alpha = alpha, burnin = burnin, m = m,
    TV_size = TV_size, nu = nu, L = L,
    Kernel = Kernel, seed = seed, n_cores = n_cores
  )

  cat("Step 1 Edge Accuracy:\n")
  print(step1$edge_accuracy)

  tv_tibble = run_node_tv_p7_parallel(
    x_list = step1$x_list,
    Edge_list = step1$Edge_list,
    k_list = step1$k_list,
    M_list_list = step1$M_list_list,
    nodes_to_test = 1:7,
    Kernel = Kernel, nu = nu, L = L,
    alpha = alpha, n_cores = n_cores
  )

  TV_results = tv_tibble |>
    group_by(node, i) |>
    summarise(
      p_min = min(p.adjust(c(Re, Im), method = "BY")),
      p_sub_used = first(p_sub),
      .groups = "drop"
    ) |>
    mutate(tv_detected = p_min < alpha)

  final = twostep_recovery_p7(step1, TV_results, alpha = alpha)

  cat("\nFinal Results:\n")
  print(final)

  return(list(
    step1 = step1,
    tv_tibble = tv_tibble,
    TV_results = TV_results,
    final = final
  ))
}

# RUN
results = run_twostep_p7(
  R = 100, alpha = 0.05, burnin = 200, m = 2^12,
  TV_size = 0.6, nu = 2, L = 1,
  Kernel = 'Kernel_Triangular', seed = 1, n_cores = 20
)
