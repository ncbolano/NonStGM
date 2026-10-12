# Running the simulations
# x = mvrnorm(2^11,Sigma = solve(A),mu=rep(0,3))
# Simulate ts
require(tidyverse)
#source('variance_functions_SSR_v5.R')
source('extractVariables.R')
source('exampleVAR.R')
source('KernelWeights.R')
source('Combined_MK_Estimation.R')
source('DFTransform.R')
source('R_Hat_Creation.R')
source('BetaFunctions.R')
source('VarianceFunctions.R')
library(doParallel)
cl<-makeCluster(16)
print(detectCores())
registerDoParallel(cl)
args=(commandArgs(TRUE))
if(length(args)==0){
  print("No arguments supplied.")
  ##supply default values
  R <- 500 # number of replications
  nu <- 2 # dimension of matrix
  p <- 3 # dimension of time series
  n <- 2^12
  M<-30
  TV_size=0.6
  L=2
  Kernel='Kernel_Triangular'
}else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}
filenames=sprintf('n_%i_M_%i_nu_%i_L_%i_TV_%0.2f_%s_extractK',n,M,nu,L,TV_size,Kernel)
################################################################
####### Stationary time series
################################################################
sim.VAR <- function(A, burnin, m) {
  p <- nrow(A)
  x <- matrix(rnorm((m + burnin) * p), ncol = p)
  x1 <- x
  for (tt in (2:(m + burnin))) {
    temp <- A %*% matrix(x1[(tt - 1), ], ncol = 1) + matrix(x[tt, ], ncol = 1)
    x1[tt, ] <- c(temp)
  }
  x2 <- x1[-c(1:burnin), ]
  return(x2)
}
################################################################
####### Nonstationary time series
################################################################
sim.tvVAR <- function(burnin, m) {
  A <- matrix(c(0.5, 0.2, 0, 0, 0.8, 0, 0, 0.3, 0.6), ncol = 3, byrow = T)
  n <- m + burnin
  p <- 3
  x <- matrix(rnorm(n * p), ncol = p)
  x1 <- x
  st <- 0.3 + TV_size * (1 + exp(0.005 * (c(1:n) - (n / 2))))**(-1)
  for (tt in (2:n)) {
    A.t <- A
    A.t[1, 1] <- st[tt]
    temp <- A.t %*% matrix(x1[tt - 1, ], ncol = 1) + matrix(x[tt, ], ncol = 1)
    x1[tt, ] <- c(temp)
  }
  x2 <- x1[-c(1:burnin), ]
  return(x2)
}
p=3
############################ Checks ############################
#######################################################
############ REAL and IMAG ############################
#######################################################
# If nu = 1, then use 8, if nu = 2 use 14
coefnum <- (p * 2 * nu + p - 1)
Kernel_function=get(Kernel)
set.seed(0)

sim_Res=foreach(i=icount(R),.combine = cbind,.inorder = FALSE,.errorhandling = "stop") %dopar% {
  n_k <- 2  # Two k values from extractK
  betaCoefAll <- array(0, c(n_k, coefnum, p))
  varbeta <- array(0, c(n_k, coefnum, 2, p))
  varbetaAll <- array(0, c(n_k, 2 * coefnum, 2 * coefnum, p))
  x <- sim.tvVAR(burnin = 20, m = n)
  JJ0 <- mvfft(x) / sqrt(nrow(x))
  n <- nrow(x)
  backward <- c(((n / 2) + 1):n)
  # we redefine JJ in this way so that we do not have to deal with the boundary close to k = 0.
  JJ <- rbind(JJ0[backward, ], JJ0[c(1:(n / 2)), ])
  delta <- 0
  #########################################
  # parameters
  J <- JJ
  # Extract k values using extractK (depends on JJ)
  ks <- extractK(J, coefnum)  # Returns two k values

  # select frequency
  # Estimating beta coefficients
  for (k in seq_along(ks)) {
    for (i_p in 1:p)
    {
      betaCoefAll[k, , i_p] <- beta(J, ks[k], nu, Kernel_function, M, a = i_p, delta = 0)
      varbetaAll[k, , , i_p] <- variance.estimator.v2(J, ks[k], nu, Kernel_function, M, a = i_p, L, delta = 0)
      varbeta[k, , 1, i_p] <- diag(varbetaAll[k, 1:coefnum, 1:coefnum, i_p])
      varbeta[k, , 2, i_p] <- diag(varbetaAll[k, 1:coefnum + coefnum, 1:coefnum + coefnum, i_p])
    }
  }
  list(betaCoefAll, varbetaAll, varbeta, ks)  # Also return ks for reference
}
betaCoefAll=sapply(sim_Res[1,],function(x) x,simplify = 'array')|>aperm(c(4,1:3))
varbetaAll=sapply(sim_Res[2,],function(x) x,simplify = 'array')|>aperm(c(5,1:4))
varbeta=sapply(sim_Res[3,],function(x) x,simplify = 'array')|>aperm(c(5,1:4))
ks_all=sapply(sim_Res[4,],function(x) x,simplify = 'array')|>t()  # R x 2 matrix of k values

r_to_loc <- function(r, c, a, nu, p) {
  if (r < 0) {
    return(p * (r + nu) + c)
  }
  if (r == 0 & c > a) {
    return(p * nu + (c - 1))
  }
  if (r == 0 & c < a) {
    return(p * nu + c)
  }
  if (r > 0) {
    return(p * nu + (p - 1) + (r - 1) * p + c)
  }
}

Test_tibble <- NULL
for (a in 1:p) {
  for (c in 1:p) {
    for (r in 0:(nu)) {
      if (!(a == c & r == 0)) {
        loc <- r_to_loc(r, c, a, nu = nu, p = p)
        # Now loop over replications since each has different ks
        for (i_rep in 1:R) {
          for (k in 1:2) {
            tmp1 <- pnorm(abs(Re(betaCoefAll[i_rep, k, loc, a])) / sqrt(abs(varbeta[i_rep, k, loc, 1, a])), lower.tail = F) * 2
            tmp2 <- pnorm(abs(Im(betaCoefAll[i_rep, k, loc, a])) / sqrt(abs(varbeta[i_rep, k, loc, 2, a])), lower.tail = F) * 2
            tmp3 <- min(p.adjust(c(tmp1, tmp2), method = "BY"))
            tmp1_vec <- c(Re(betaCoefAll[i_rep, k, loc, a]), Im(betaCoefAll[i_rep, k, loc, a]))
            matrix_tmp <- varbetaAll[i_rep, k, c(loc, loc + coefnum), c(loc, loc + coefnum), a]
            matrix_tmp[1, 2] <- 0
            matrix_tmp[2, 1] <- 0
            matrix_tmp <- abs(matrix_tmp)
            tmp4 <- pchisq(tmp1_vec %*% solve(matrix_tmp, b = tmp1_vec), lower.tail = F, df = 2)

            Test_tibble <- tibble(Re = tmp1, Im = tmp2, padjust = tmp3, Chisq = tmp4,
                                  a = a, c = c, r = r, k = ks_all[i_rep, k], k_idx = k, i = i_rep) |>
              rbind(Test_tibble)
          }
        }
      }
    }
  }
}

Test_tibble|>mutate(n=n,M=M,nu=nu,TV_size=TV_size,Kernel=Kernel)

Graph_tibble <- Test_tibble |>
  filter(a <= c) |>  # Only consider a <= c
  group_by(a, c, i) |>
  summarise(
    p_min = min(p.adjust(c(Re, Im), method = "BY")),
    .groups = "drop"
  ) |>
  select(i, a, c, p_min)

alpha = 0.05

# Option 3: Binary adjacency indicator per replication
Recovered_binary <- Graph_tibble |>
  mutate(edge = as.integer(p_min < alpha)) |>
  select(i, a, c, edge) |>
  pivot_wider(
    names_from = c(a, c),
    values_from = edge,
    names_prefix = "edge_"
  )



true_graph <- c(
  edge_1_1 = 1,  # A[1,1] = TV
  edge_1_2 = 1,  # A[1,2] = 0.2
  edge_1_3 = 0,
  edge_2_2 = 0,
  edge_2_3 = 1,  # A[2,3] = 0.3
  edge_3_3 = 0
)

n_edges <- length(true_graph)
n_true <- sum(true_graph)
n_null <- sum(true_graph == 0)

# Extract just the edge columns (remove i)
recovered_matrix <- Recovered_binary |>
  select(-i) |>
  as.matrix()

# Metrics per replication
Metrics <- tibble(
  i = Recovered_binary$i,

  TP = rowSums(sweep(recovered_matrix, 2, true_graph, `*`)),

  FP = rowSums(sweep(recovered_matrix, 2, 1 - true_graph, `*`)),

  FN = rowSums(sweep(1 - recovered_matrix, 2, true_graph, `*`)),


  TN = rowSums(sweep(1 - recovered_matrix, 2, 1 - true_graph, `*`))
) |>
  mutate(
    perfect_recovery = (TP == n_true) & (FP == 0),

    power = TP / n_true,

    FDR = ifelse(TP + FP == 0, 0, FP / (TP + FP)),

    specificity = TN / n_null
  )

# Overall summary
Metrics |>
  summarise(
    perfect_recovery_rate = mean(perfect_recovery),
    avg_power = mean(power),
    avg_FDR = mean(FDR),
    avg_specificity = mean(specificity),
    avg_TP = mean(TP),
    avg_FP = mean(FP)
  )

