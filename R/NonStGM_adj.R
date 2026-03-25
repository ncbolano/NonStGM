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

alpha = .05
tib = NonStGM_adj(x, Kernel = 'Kernel_Triangular', nu = 2, L = 1, alpha = alpha)
graph_recovery(tib, alpha = alpha)


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
