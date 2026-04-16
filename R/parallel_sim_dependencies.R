
# ============================================================
# p=3 VAR SIMULATION FUNCTION
# ============================================================
source('extractVariables.R')
source('exampleVAR.R')
source('KernelWeights.R')
source('Combined_MK_Estimation.R')
source('DFTransform.R')
source('R_Hat_Creation.R')
source('BetaFunctions.R')
source('VarianceFunctions.R')
source('NonStGM_adj.R')

library(foreach)
library(doParallel)
library(parallel)
library(tidyverse)

p = 3

sim.tvVAR_p3 = function(burnin, m, TV_size) {
  A = matrix(c(
    0.5, 0.2, 0,
    0,   0.8, 0,
    0,   0.3, 0.6
  ), ncol = 3, byrow = TRUE)

  n = m + burnin
  p = 3
  x = matrix(rnorm(n * p), ncol = p)
  x1 = x
  st = 0.3 + TV_size * (1 + exp(0.005 * (c(1:n) - (n / 2))))^(-1)

  for (tt in 2:n) {
    A.t = A
    A.t[1, 1] = st[tt]
    temp = A.t %*% matrix(x1[tt - 1, ], ncol = 1) + matrix(x[tt, ], ncol = 1)
    x1[tt, ] = c(temp)
  }

  x2 = x1[-c(1:burnin), ]
  return(x2)
}

check_tvVAR_stability_p3 = function(burnin, m, TV_size) {
  A = matrix(c(
    0.5, 0.2, 0,
    0,   0.8, 0,
    0,   0.3, 0.6
  ), ncol = 3, byrow = TRUE)

  n = burnin + m
  st = 0.3 + TV_size * (1 + exp(0.005 * (1:n - n/2)))^(-1)
  max_mod = sapply(st, function(s) { A[1, 1] = s; max(Mod(eigen(A)$values)) })
  cat("A[1,1] range:", round(range(st), 3), "| Max eigenvalue modulus:", round(max(max_mod), 4),
      "| Stable:", max(max_mod) < 1, "\n")
  return(max(max_mod) < 1)
}

# ============================================================
# GRAPH RECOVERY FUNCTION FOR p=3
# ============================================================
# True edges:
#   Time-varying: (1,1)
#   Off-diagonal: (1,2), (2,3)
#   Total: 3 true edges, 2 off-diagonal

graph_recovery_p3 = function(Test_tibble, alpha = 0.05) {

  Test_tibble = Test_tibble |>
    filter(a <= c) |>
    filter((a != c & r == 0) | (a == c & r != 0))

  Graph_tibble = Test_tibble |>
    group_by(a, c, i) |>
    summarise(
      p_min = min(p.adjust(c(Re, Im), method = "BY")),
      .groups = "drop"
    )

  graph_recovery_results = Graph_tibble |>
    group_by(i) |>
    summarise(
      # True edges
      detect_11_tv = any(a == 1 & c == 1 & p_min < alpha),
      detect_12 = any(a == 1 & c == 2 & p_min < alpha),
      detect_23 = any(a == 2 & c == 3 & p_min < alpha),
      # False edges (null)
      detect_13 = any(a == 1 & c == 3 & p_min < alpha),
      detect_22 = any(a == 2 & c == 2 & p_min < alpha),
      detect_33 = any(a == 3 & c == 3 & p_min < alpha),
      # Edge counts
      n_offdiag_edges = sum(a != c & p_min < alpha),
      n_diag_edges = sum(a == c & p_min < alpha),
      .groups = "drop"
    )

  accuracy_summary = graph_recovery_results |>
    mutate(
      all_offdiag_detected = detect_12 & detect_23,
      all_true_detected = all_offdiag_detected & detect_11_tv,
      perfect = all_offdiag_detected & (n_offdiag_edges == 2),
      sensitivity = (detect_11_tv + detect_12 + detect_23) / 3,
      false_positives = detect_13 + detect_22 + detect_33
    )

  overall_accuracy = tibble(
    Method = "Global BY",
    Perfect_Recovery = mean(accuracy_summary$perfect),
    All_True_Detected = mean(accuracy_summary$all_true_detected),
    # True edges
    Detect_11_TV = mean(accuracy_summary$detect_11_tv),
    Detect_12 = mean(accuracy_summary$detect_12),
    Detect_23 = mean(accuracy_summary$detect_23),
    # False edges
    False_13 = mean(accuracy_summary$detect_13),
    False_22 = mean(accuracy_summary$detect_22),
    False_33 = mean(accuracy_summary$detect_33),
    # Summary
    Mean_Sensitivity = mean(accuracy_summary$sensitivity),
    Mean_N_Offdiag_Edges = mean(accuracy_summary$n_offdiag_edges),
    Mean_False_Positives = mean(accuracy_summary$false_positives)
  )

  return(overall_accuracy)
}


p = 5

sim.tvVAR_p5 = function(burnin, m, TV_size) {
  A = matrix(c(
    0.5,  0,    0,    0,    0,
    0.3,  0.6,  0,    0,    0,
    0,    0.3,  0.5,  0,    0,
    0,    0,    0.3, 0.55, 0,
    0.3, 0,    0,    0.3,  0.5
  ), ncol = 5, byrow = TRUE)
  n = m + burnin
  p = 5
  x = matrix(rnorm(n * p), ncol = p)
  x1 = x
  st = 0.3 + TV_size * (1 + exp(0.005 * (c(1:n) - (n / 2))))^(-1)
  for (tt in 2:n) {
    A.t = A
    A.t[1, 1] = st[tt]
    temp = A.t %*% matrix(x1[tt - 1, ], ncol = 1) + matrix(x[tt, ], ncol = 1)
    x1[tt, ] = c(temp)
  }
  x2 = x1[-c(1:burnin), ]
  return(x2)
}

check_tvVAR_stability = function(burnin, m, TV_size) {
  A = matrix(c(
    0.5,  0,    0,    0,    0,
    0.3,  0.6,  0,    0,    0,
    0,    0.3,  0.5,  0,    0,
    0,    0,    0.3, 0.55, 0,
    0.3, 0,    0,    0.3,  0.5
  ), ncol = 5, byrow = TRUE)
  n = burnin + m
  st = 0.3 + TV_size * (1 + exp(0.005 * (1:n - n/2)))^(-1)
  max_mod = sapply(st, function(s) { A[1,1] = s; max(Mod(eigen(A)$values)) })
  cat("A[1,1] range:", round(range(st), 3), "| Max eigenvalue modulus:", round(max(max_mod), 4),
      "| Stable:", max(max_mod) < 1, "\n")
  return(max(max_mod) < 1)
}

graph_recovery_p5 = function(Test_tibble, alpha = 0.05) {
  Test_tibble = Test_tibble |>
    filter(a <= c) |>
    filter((a != c & r == 0) | (a == c & r != 0))

  Graph_tibble = Test_tibble |>
    group_by(a, c, i) |>
    summarise(
      p_min = min(p.adjust(c(Re, Im), method = "BY")),
      .groups = "drop"
    )

  graph_recovery_results = Graph_tibble |>
    group_by(i) |>
    summarise(
      detect_11_tv = any(a == 1 & c == 1 & p_min < alpha),
      detect_12 = any(a == 1 & c == 2 & p_min < alpha),
      detect_23 = any(a == 2 & c == 3 & p_min < alpha),
      detect_34 = any(a == 3 & c == 4 & p_min < alpha),
      detect_45 = any(a == 4 & c == 5 & p_min < alpha),
      detect_15 = any(a == 1 & c == 5 & p_min < alpha),
      n_offdiag_edges = sum(a != c & p_min < alpha),
      .groups = "drop"
    )

  accuracy_summary = graph_recovery_results |>
    mutate(
      all_offdiag_detected = detect_12 & detect_23 & detect_34 &
        detect_45 & detect_15,
      all_true_detected = all_offdiag_detected & detect_11_tv,
      perfect = all_offdiag_detected & (n_offdiag_edges == 5),
      sensitivity = (detect_11_tv + detect_12 + detect_23 +
                       detect_34 + detect_45 + detect_15) / 6,
      false_positives = pmax(0, n_offdiag_edges - 5)
    )

  overall_accuracy = tibble(
    Method = "Global BY",
    Perfect_Recovery = mean(accuracy_summary$perfect),
    All_True_Detected = mean(accuracy_summary$all_true_detected),
    Detect_11_TV = mean(accuracy_summary$detect_11_tv),
    Detect_12 = mean(accuracy_summary$detect_12),
    Detect_23 = mean(accuracy_summary$detect_23),
    Detect_34 = mean(accuracy_summary$detect_34),
    Detect_45 = mean(accuracy_summary$detect_45),
    Detect_15 = mean(accuracy_summary$detect_15),
    Mean_Sensitivity = mean(accuracy_summary$sensitivity),
    Mean_N_Edges = mean(accuracy_summary$n_offdiag_edges),
    Mean_False_Positives = mean(accuracy_summary$false_positives)
  )
  return(overall_accuracy)
}

# ============================================================
# SIMULATION RUNNER FOR p=5
# ============================================================
run_NonStGM_simulation_p5 = function(R = 100, alpha = 0.05, burnin = 20, m = 4096,
                                     TV_size = 0.6, nu = 2, L = 1,
                                     Kernel = 'Kernel_Triangular', seed = 0) {
  All_Test_tibble = NULL
  set.seed(seed)
  for (iter in 1:R) {
    x = sim.tvVAR_p5(burnin, m, TV_size)
    tib = NonStGM_adj(x, Kernel = Kernel, nu = nu, L = L, alpha = alpha)
    tib = tib |> mutate(i = iter)
    All_Test_tibble = rbind(All_Test_tibble, tib)
  }
  results = graph_recovery_p5(All_Test_tibble, alpha = alpha)
  return(results)
}

# results_p5 = run_NonStGM_simulation_p5(
#   R = 100,
#   alpha = 0.05,
#   burnin = 200,
#   m = 2^12,
#   TV_size = 0.6,
#   nu = 2,
#   L = 1,
#   Kernel = 'Kernel_Triangular',
#   seed = 1
# )
# print(results_p5)

# ============================================================
# p=7 VAR SIMULATION FUNCTION
# ============================================================
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
  st = .3 + TV_size * (1 + exp(0.005 * (c(1:n) - (n / 2))))^(-1) # CHANGED
  for (tt in 2:n) {
    A.t = A
    A.t[1, 1] = st[tt]
    temp = A.t %*% matrix(x1[tt - 1, ], ncol = 1) + matrix(x[tt, ], ncol = 1)
    x1[tt, ] = c(temp)
  }
  x2 = x1[-c(1:burnin), ]
  return(x2)
}

check_tvVAR_stability_p7 = function(burnin, m, TV_size) {
  A = matrix(c(
    0.5,  0,    0,    0,    0,    0,    0,
    0.3,  0.5,  0,    0,    0,    0,    0,
    0,    0.3,  0.5,  0,    0,    0,    0,
    0,    0,    0.3,  0.5,  0,    0,    0,
    0,    0,    0,    0.3,  0.5,  0,    0,
    0,    0,    0,    0,    0.3,  0.5,  0,
    0.3,  0,    0,    0,    0,    0.3,  0.5
  ), ncol = 7, byrow = TRUE)
  n = burnin + m
  st = .3 + TV_size * (1 + exp(0.005 * (1:n - n/2)))^(-1) #CHANGED
  max_mod = sapply(st, function(s) { A[1, 1] = s; max(Mod(eigen(A)$values)) })
  cat("A[1,1] range:", round(range(st), 3), "| Max eigenvalue modulus:", round(max(max_mod), 4),
      "| Stable:", max(max_mod) < 1, "\n")
  return(max(max_mod) < 1)
}


graph_recovery_p7 = function(Test_tibble, alpha = 0.05) {
  Test_tibble = Test_tibble |>
    filter(a <= c) |>
    filter((a != c & r == 0) | (a == c & r != 0))

  Graph_tibble = Test_tibble |>
    group_by(a, c, i) |>
    summarise(
      p_min = min(p.adjust(c(Re, Im), method = "BY")),
      .groups = "drop"
    )

  graph_recovery_results = Graph_tibble |>
    group_by(i) |>
    summarise(
      detect_11_tv = any(a == 1 & c == 1 & p_min < alpha),
      detect_12 = any(a == 1 & c == 2 & p_min < alpha),
      detect_23 = any(a == 2 & c == 3 & p_min < alpha),
      detect_34 = any(a == 3 & c == 4 & p_min < alpha),
      detect_45 = any(a == 4 & c == 5 & p_min < alpha),
      detect_56 = any(a == 5 & c == 6 & p_min < alpha),
      detect_67 = any(a == 6 & c == 7 & p_min < alpha),
      detect_17 = any(a == 1 & c == 7 & p_min < alpha),
      n_offdiag_edges = sum(a != c & p_min < alpha),
      .groups = "drop"
    )

  accuracy_summary = graph_recovery_results |>
    mutate(
      all_offdiag_detected = detect_12 & detect_23 & detect_34 &
        detect_45 & detect_56 & detect_67 & detect_17,
      all_true_detected = all_offdiag_detected & detect_11_tv,
      perfect = all_offdiag_detected & (n_offdiag_edges == 7),
      sensitivity = (detect_11_tv + detect_12 + detect_23 + detect_34 +
                       detect_45 + detect_56 + detect_67 + detect_17) / 8,
      false_positives = pmax(0, n_offdiag_edges - 7)
    )

  overall_accuracy = tibble(
    Method = "Global BY",
    Perfect_Recovery = mean(accuracy_summary$perfect),
    All_True_Detected = mean(accuracy_summary$all_true_detected),
    Detect_11_TV = mean(accuracy_summary$detect_11_tv),
    Detect_12 = mean(accuracy_summary$detect_12),
    Detect_23 = mean(accuracy_summary$detect_23),
    Detect_34 = mean(accuracy_summary$detect_34),
    Detect_45 = mean(accuracy_summary$detect_45),
    Detect_56 = mean(accuracy_summary$detect_56),
    Detect_67 = mean(accuracy_summary$detect_67),
    Detect_17 = mean(accuracy_summary$detect_17),
    Mean_Sensitivity = mean(accuracy_summary$sensitivity),
    Mean_N_Edges = mean(accuracy_summary$n_offdiag_edges),
    Mean_False_Positives = mean(accuracy_summary$false_positives)
  )
  return(overall_accuracy)
}

# ============================================================
# SIMULATION RUNNER FOR p=7
# ============================================================
run_NonStGM_simulation_p7 = function(R = 100, alpha = 0.05, burnin = 20, m = 4096,
                                     TV_size = 0.6, nu = 2, L = 1,
                                     Kernel = 'Kernel_Triangular', seed = 0) {
  All_Test_tibble = NULL
  set.seed(seed)
  for (iter in 1:R) {
    x = sim.tvVAR_p7(burnin, m, TV_size)
    tib = NonStGM_adj(x, Kernel = Kernel, nu = nu, L = L, alpha = alpha)
    tib = tib |> mutate(i = iter)
    All_Test_tibble = rbind(All_Test_tibble, tib)
  }
  results = graph_recovery_p7(All_Test_tibble, alpha = alpha)
  return(list(results, All_Test_tibble))
}

# results_p7 = run_NonStGM_simulation_p7(
#   R = 100,
#   alpha = 0.05,
#   burnin = 20,
#   m = 2^12,
#   TV_size = 0.6,
#   nu = 2,
#   L = 1,
#   Kernel = 'Kernel_Triangular',
#   seed = 0
# )
# print(results_p7[[1]])
# hi = results_p7[[2]]
