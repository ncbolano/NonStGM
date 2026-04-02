# Running the simulations
# x = mvrnorm(2^11,Sigma = solve(A),mu=rep(0,3))
# Simulate ts
require(tidyverse)
source('variance_functions_SSR_v5.R')


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
  L=1
  ks <- c((n / 2 + M), (n / 2 + M + 2 * M))
  Kernel='Kernel_Triangular'

}else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}


filenames=sprintf('n_%i_M_%i_nu_%i_L_%i_TV_%0.2f_%s',n,M,nu,L,TV_size,Kernel)

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
# Ms=array(0,c(R,n_k,p))
# M_grid <- c(20 + 0:5 * 30)

set.seed(0)

sim_Res=foreach(i=icount(R),.combine = cbind,.inorder = FALSE,.errorhandling = "stop") %dopar% {

  n_k <- length(ks)
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


  # select frequency



  # Estimating beta coefficients

  for (k in seq_along(ks)) {
    for (i_p in 1:p)
    {
      # M <- M_grid[which.min(beta_cv2(J, ks[k], nu,Kernel_Triangular, M_grid = M_grid, a = i_p, delta = 0))]
      # Ms[i,k,i_p]<- M
      betaCoefAll[k, , i_p] <- beta(J, ks[k], nu,Kernel_function, M, a = i_p, delta = 0)

      varbetaAll[k, , , i_p] <- variance.estimator.v2(J, ks[k], nu,Kernel_function, M, a = i_p, L, delta = 0)

      varbeta[k, , 1, i_p] <- diag(varbetaAll[k, 1:coefnum, 1:coefnum, i_p])
      varbeta[k, , 2, i_p] <- diag(varbetaAll[k, 1:coefnum + coefnum, 1:coefnum + coefnum, i_p])
    }
  }
  list(betaCoefAll,varbetaAll,varbeta)
}

betaCoefAll=sapply(sim_Res[1,],function(x) x,simplify = 'array')|>aperm(c(4,1:3))
varbetaAll=sapply(sim_Res[2,],function(x) x,simplify = 'array')|>aperm(c(5,1:4))
varbeta=sapply(sim_Res[3,],function(x) x,simplify = 'array')|>aperm(c(5,1:4))
# save(betaCoefAll,varbeta,varbetaAll,ks,file=sprintf('Simu_Data_%s.Rdata',filenames))

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
      for (k in seq_along(ks))
      {
        if (!(a == c & r == 0)) {
          loc <- r_to_loc(r, c, a, nu = nu, p = p)
          tmp1 <- pnorm(abs(Re(betaCoefAll[, k, loc, a])) / sqrt(abs(varbeta[, k, loc, 1, a])), lower.tail = F) * 2
          tmp2 <- pnorm(abs(Im(betaCoefAll[, k, loc, a])) / sqrt(abs(varbeta[, k, loc, 2, a])), lower.tail = F) * 2
          tmp3 <- unlist(purrr::map2(tmp1, tmp2, function(x, y) min(p.adjust(c(x, y), method = "BY"))))
          tmp4 <- sapply(1:R, function(i) {
            tmp1 <- c(Re(betaCoefAll[i, k, loc, a]), Im(betaCoefAll[i, k, loc, a]))

            matrix_tmp <- varbetaAll[i, k, c(loc, loc + coefnum), c(loc, loc + coefnum), a]
            matrix_tmp[1, 2] <- 0
            matrix_tmp[2, 1] <- 0
            matrix_tmp <- abs(matrix_tmp)
            return(pchisq(tmp1 %*% solve(matrix_tmp, b = tmp1), lower.tail = F, df = 2))
          })
          Test_tibble <- cbind(Re = tmp1, Im = tmp2, padjust = tmp3, Chisq = tmp4) |>
            as_tibble() |>
            mutate(a = a, c = c, r = r, k = ks[k], i = 1:R) |>
            rbind(Test_tibble)
        }
      }
    }
  }
}

Test_tibble|>mutate(n=n,M=M,nu=nu,TV_size=TV_size,Kernel=Kernel)|>saveRDS(file=sprintf('Simu_tbl_%s.RData',filenames))
#
# alpha <- 0.05
# erg1 <- Test_tibble |>
#   group_by(a, c, r, k) |>
#   summarise(across(-i, ~ mean(.x < alpha))) |>
#   arrange(k, a, c, r)
#
#
 erg2 <- Test_tibble |>
   mutate(tv = r > 0) |> # Creates new columns of test_tibble for all r > 0
   group_by(a, c, tv, i) |> # Groups by our four variables
   summarise(
     padjust2 = min(p.adjust(c(Re, Im), method = "BY")),
     padjust3 = min(p.adjust(c(Re, Im), method = "none"))
   ) |>
   group_by(a, c, tv) |>
   summarise(across(everything(), ~ mean(.x < alpha))) |>
   arrange(a, c, tv)

#
# erg1|>filter(k==ks[1],a<=c)|>left_join(erg2,join_by(a,c,r))|>print(n=300)
#
#
# erg1|>filter(k==ks[2],a<=c)|>left_join(erg2,join_by(a,c,r))|>print(n=300)

?p.adjust

# a=2
# c=2
# r=2
# qqnorm(y=(Re(betaCoefAll[, r_to_loc(r, c, a,nu = nu,p = p), a])) / sqrt(abs(varbeta[, r_to_loc(r, c, a,nu = nu,p = p), 1, a])))
#
#
# (2 * M + 1) * apply(Re(betaCoefAll[,,1]), 2, var)
# (2 * M + 1) * apply(varbeta[,,1,1], 2, mean)
# (2 * M + 1) * apply(Im(betaCoefAll[,,1]), 2, var)
# (2 * M + 1) * apply(varbeta[,,2,1], 2, mean)
#
#
#..
# a=1
# c=1
# r=2
# tmp1=(2 * M + 1) *var(cbind(Re(betaCoefAll[,,1]),Im(betaCoefAll[,,a])))
# tmp2=(2 * M + 1) *apply(varbetaAll[,,,a],2:3,mean)
#
# loc <- r_to_loc(r, c, a,nu = nu,p = p)
# tmp1[c(loc, c(loc) + coefnum),c(loc, c(loc) + coefnum)]
# tmp2[c(loc, c(loc) + coefnum),c(loc, c(loc) + coefnum)]
# hist((2*M+1)*varbetaAll[,loc,loc,a],breaks=25)
# abline(v=(2*M+1)*mean(varbetaAll[,loc,loc,a]),col='red')

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
