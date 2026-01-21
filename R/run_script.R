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

# sim function
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
list(betaCoefAll, varbetaAll, varbeta)

betaCoefAll = sapply(sim_Res[1, ], function(x) x, simplify = 'array') |> aperm(c(4, 1:3))
varbetaAll = sapply(sim_Res[2, ], function(x) x, simplify = 'array') |> aperm(c(5, 1:4))
varbeta = sapply(sim_Res[3, ], function(x) x, simplify = 'array') |> aperm(c(5, 1:4))

Test_tibble = NULL
for (a in 1:p) {
  for (c in 1:p) {
    for (r in 0:(nu)) {
      for (j in seq_along(k)) {
        if (!(a == c & r == 0)) {
          loc = r_to_loc(r, c, a, nu = nu, p = p)
          tmp1 = pnorm(abs(Re(betaCoefAll[, j, loc, a])) / sqrt(abs(varbeta[, j, loc, 1, a])), lower.tail = F) * 2
          tmp2 = pnorm(abs(Im(betaCoefAll[, j, loc, a])) / sqrt(abs(varbeta[, j, loc, 2, a])), lower.tail = F) * 2
          tmp3 = unlist(purrr::map2(tmp1, tmp2, function(x, y) min(p.adjust(c(x, y), method = "BY"))))
          tmp4 = sapply(1:R, function(i) {
            tmp1_vec = c(Re(betaCoefAll[i, j, loc, a]), Im(betaCoefAll[i, j, loc, a]))
            matrix_tmp = varbetaAll[i, j, c(loc, loc + coefnum), c(loc, loc + coefnum), a]
            matrix_tmp[1, 2] = 0  # Fixed indexing
            matrix_tmp[2, 1] = 0  # Fixed indexing
            matrix_tmp = abs(matrix_tmp)
            return(pchisq(tmp1_vec %*% solve(matrix_tmp, b = tmp1_vec), lower.tail = F, df = 2))
          })
          Test_tibble = cbind(Re = tmp1, Im = tmp2, padjust = tmp3, Chisq = tmp4) |>
            as_tibble() |>
            mutate(a = a, c = c, r = r, k1 = k[j], i = 1:R) |>
            rbind(Test_tibble)
        }
      }
    }
  }
}

