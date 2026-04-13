NonStGM_adj = function(x, Kernel = 'Kernel_Triangular', nu = 2, L = 1, alpha = alpha) {

  Kernel_function = get(Kernel)

  # Load DFT, dimension variables, coefficient number
  variable_list = extractVariables(x)
  J = variable_list[[1]] ; p = as.numeric(variable_list[2])
  coefnum = (2 * p * nu) + p - 1

  # Gathering frequencies (k) + local smoothing values (m_i) for each frequency
  k = extractK(J, coefnum)
  n = nrow(x)
  M_grid = unique(round(seq(max(coefnum , n^(1/5)), max(coefnum , n^(1/2))  , length.out = 50)))
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
      for (r in 0:(nu)) {
        for (j in seq_along(k)) {
          if (!(a == c & r == 0)) {
            loc = r_to_loc(r, c, a, nu = nu, p = p)

            tmp1 = pnorm(abs(Re(betaCoefAll[j, loc, a])) / sqrt(abs(varbeta[j, loc, 1, a])), lower.tail = F) * 2
            tmp2 = pnorm(abs(Im(betaCoefAll[j, loc, a])) / sqrt(abs(varbeta[j, loc, 2, a])), lower.tail = F) * 2
            tmp3 = min(p.adjust(c(tmp1, tmp2), method = "BY"))

            tmp1_vec = c(Re(betaCoefAll[j, loc, a]), Im(betaCoefAll[j, loc, a]))
            matrix_tmp = varbetaAll[j, c(loc, loc + coefnum), c(loc, loc + coefnum), a]
            matrix_tmp[1, 2] = 0
            matrix_tmp[2, 1] = 0
            matrix_tmp = abs(matrix_tmp)
            tmp4 = pchisq(tmp1_vec %*% solve(matrix_tmp, b = tmp1_vec), lower.tail = F, df = 2)

            Test_tibble = tibble(
              Re = tmp1,
              Im = tmp2,
              padjust = tmp3,
              Chisq = tmp4,
              a = a,
              c = c,
              r = r,
              k = k[j],
              i = 1
            ) |>
              rbind(Test_tibble)
          }
        }
      }
    }
  }
  Test_tibble = Test_tibble |>
    filter(a <= c)
  return(Test_tibble)
}

#alpha = .05
#tib = NonStGM_adj(x, Kernel = 'Kernel_Triangular', nu = 2, L = 1, alpha = alpha)

# ============================================================================
# GRAPH RECOVERY FUNCTION
# ============================================================================

# graph_recovery = function(Test_tibble, alpha = 0.05) {
#
#
#   # Filter to lower triangular
#   #Test_tibble = Test_tibble |>
#     #filter(a <= c)
#
#   # Global BY adjustment
#   Test_tibble_global = Test_tibble |>
#     group_by(i) |>
#     group_modify(\(df, key) {
#       all_pvals = c(df$Re, df$Im)
#       all_adj = p.adjust(all_pvals, method = "BY")
#       n_tests = nrow(df)
#       df$Re_global_BY = all_adj[1:n_tests]
#       df$Im_global_BY = all_adj[(n_tests + 1):(2 * n_tests)]
#       df$global_BY_min = pmin(df$Re_global_BY, df$Im_global_BY)
#       df
#     }) |>
#     ungroup()
#
#   # Two-step BY adjustment
#   # Step 1: Within each (a, c, r, k, i) group
#   Test_tibble_twostep = Test_tibble |>
#     group_by(a, c, r, k, i) |>
#     group_modify(\(df, key) {
#       all_pvals = c(df$Re, df$Im)
#       all_adj = p.adjust(all_pvals, method = "BY")
#       n_tests = nrow(df)
#       df$Re_step1_BY = all_adj[1:n_tests]
#       df$Im_step1_BY = all_adj[(n_tests + 1):(2 * n_tests)]
#       df$step1_min = pmin(df$Re_step1_BY, df$Im_step1_BY)
#       df
#     }) |>
#     ungroup()
#
#   # Step 2: Apply BY correction across all step1_min values
#   Test_tibble_twostep = Test_tibble_twostep |>
#     mutate(twostep_BY = p.adjust(step1_min, method = "BY"))
#
#   # Step 3: Merge results
#   Test_tibble_final = Test_tibble |>
#     left_join(
#       Test_tibble_global |> select(a, c, r, k, i, global_BY_min),
#       by = c("a", "c", "r", "k", "i")
#     ) |>
#     left_join(
#       Test_tibble_twostep |> select(a, c, r, k, i, twostep_BY),
#       by = c("a", "c", "r", "k", "i")
#     )
#
#   # Graph recovery metrics
#   graph_recovery_results = Test_tibble_final |>
#     group_by(i) |>
#     summarise(
#       # For Global BY method
#       global_detect_12 = any(a == 1 & c == 2 & global_BY_min < alpha) |
#         any(a == 2 & c == 1 & global_BY_min < alpha),
#       global_detect_23 = any(a == 2 & c == 3 & global_BY_min < alpha) |
#         any(a == 3 & c == 2 & global_BY_min < alpha),
#       global_detect_11_tv = any(a == 1 & c == 1 & r > 0 & global_BY_min < alpha),
#       # Count number of significant edges (excluding diagonal)
#       global_n_edges = nrow(distinct(data.frame(
#         a2 = pmin(a[a != c & global_BY_min < alpha],
#                   c[a != c & global_BY_min < alpha]),
#         c2 = pmax(a[a != c & global_BY_min < alpha],
#                   c[a != c & global_BY_min < alpha])
#       ))),
#       # For Two-Step BY method
#       twostep_detect_12 = any(a == 1 & c == 2 & twostep_BY < alpha) |
#         any(a == 2 & c == 1 & twostep_BY < alpha),
#       twostep_detect_23 = any(a == 2 & c == 3 & twostep_BY < alpha) |
#         any(a == 3 & c == 2 & twostep_BY < alpha),
#       twostep_detect_11_tv = any(a == 1 & c == 1 & r > 0 & twostep_BY < alpha),
#       # Count number of significant edges for two-step
#       twostep_n_edges = nrow(distinct(data.frame(
#         a2 = pmin(a[a != c & twostep_BY < alpha],
#                   c[a != c & twostep_BY < alpha]),
#         c2 = pmax(a[a != c & twostep_BY < alpha],
#                   c[a != c & twostep_BY < alpha])
#       ))),
#       .groups = "drop"
#     )
#
#   # Accuracy summary
#   accuracy_summary = graph_recovery_results |>
#     mutate(
#       # Global BY: Perfect recovery
#       global_all_true_detected = global_detect_12 & global_detect_23 & global_detect_11_tv,
#       global_perfect = global_all_true_detected & (global_n_edges == 2),
#       # Two-Step BY: Perfect recovery
#       twostep_all_true_detected = twostep_detect_12 & twostep_detect_23 & twostep_detect_11_tv,
#       twostep_perfect = twostep_all_true_detected & (twostep_n_edges == 2),
#       # Partial credit metrics
#       global_sensitivity = (global_detect_12 + global_detect_23 + global_detect_11_tv) / 3,
#       twostep_sensitivity = (twostep_detect_12 + twostep_detect_23 + twostep_detect_11_tv) / 3
#     )
#
#   # Overall accuracy metrics
#   overall_accuracy = tibble(
#     Method = c("Global BY", "Two-Step BY"),
#     # Perfect graph recovery rate
#     Perfect_Recovery = c(
#       mean(accuracy_summary$global_perfect),
#       mean(accuracy_summary$twostep_perfect)
#     ),
#     # All true edges detected (but can still have false positives)
#     All_True_Detected = c(
#       mean(accuracy_summary$global_all_true_detected),
#       mean(accuracy_summary$twostep_all_true_detected)
#     ),
#     # Edge detection rates by pair
#     Detect_12 = c(
#       mean(accuracy_summary$global_detect_12),
#       mean(accuracy_summary$twostep_detect_12)
#     ),
#     Detect_23 = c(
#       mean(accuracy_summary$global_detect_23),
#       mean(accuracy_summary$twostep_detect_23)
#     ),
#     Detect_11_TV = c(
#       mean(accuracy_summary$global_detect_11_tv),
#       mean(accuracy_summary$twostep_detect_11_tv)
#     ),
#     # Sensitivity (proportion of true edges detected)
#     Mean_Sensitivity = c(
#       mean(accuracy_summary$global_sensitivity),
#       mean(accuracy_summary$twostep_sensitivity)
#     ),
#     # Mean number of detected edges
#     Mean_N_Edges = c(
#       mean(accuracy_summary$global_n_edges),
#       mean(accuracy_summary$twostep_n_edges)
#     ),
#     # False positive rate (extra edges attributed)
#     Mean_False_Positives = c(
#       mean(pmax(0, accuracy_summary$global_n_edges - 2)),
#       mean(pmax(0, accuracy_summary$twostep_n_edges - 2))
#     )
#   )
#
#   return(overall_accuracy)
# }

# ============================================================
# WRAPPED FUNCTION FOR R REPLICATIONS
# ============================================================

run_NonStGM_simulation <- function(R = 100, alpha = 0.05, burnin = 20, m = 4096,
                                   TV_size = 0.6, nu = 2, L = 1,
                                   Kernel = 'Kernel_Triangular', seed = 0) {

  # Collect all Test_tibbles across replications
  All_Test_tibble <- NULL

  set.seed(seed)

  for (iter in 1:R) {

    # Generate data
    x <- sim.tvVAR(burnin, m, TV_size)

    # Get Test_tibble from NonStGM_adj
    tib <- NonStGM_adj(x, Kernel = Kernel, nu = nu, L = L, alpha = alpha)

    # Update iteration index to reflect replication number
    tib <- tib |>
      mutate(i = iter)

    # Combine with all replications
    All_Test_tibble <- rbind(All_Test_tibble, tib)

    if (iter %% 20 == 0) cat("Completed iteration:", iter, "/", R, "\n")
  }

  # Apply graph_recovery to combined tibble
  results <- graph_recovery(All_Test_tibble, alpha = alpha)

  return(results)
}

# Run it
results <- run_NonStGM_simulation(R = 100, alpha = 0.05, burnin = 20, m = 4096,
                                  TV_size = 0.6, seed = 0)
print(results)
