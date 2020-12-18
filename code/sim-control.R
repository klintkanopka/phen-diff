library(tidyverse)
library(parallel)

sim.data <- function(N=1e5, rho=0.5, rho_g=0.3, p_allele=0.3,
                     alpha=-1, beta=0.04, gamma_1=sqrt(0.35),   gamma_2=sqrt(0.649),
                     continuous=FALSE){

  sim.genes <- function(N=N, rho=rho, rho_g=rho_g, p_allele=p_allele){
    g_m <- rbinom(n=N, size=2, prob=p_allele)
    g_p <- rbinom(n=N, size=2, prob=p_allele)

    out <- data.frame(g_m = g_m,
                      g_p = g_p,
                      g0 = rbinom(N, 1, g_m/2) + rbinom(N, 1, g_m/2),
                      g1 = rbinom(N, 1, g_m/2) + rbinom(N, 1, g_m/2),
                      g2 = rbinom(N, 1, g_m/2) + rbinom(N, 1, g_m/2),
                      eps_fam = rnorm(N))
    return(out)
  }

  sigmoid <- function(z){
    out <- 1 / (1 + exp(-z))
  }

  a <- rep(alpha, N)
  b <- rep(beta, N)

  d <- sim.genes(N=N, rho=rho, rho_g=rho_g, p_allele=p_allele)

  if (continuous){
    d$Y0 <- a + b*d$g0 + gamma_1*d$eps_fam + gamma_2*rnorm(N)
    d$Y1 <- a + b*d$g1 + gamma_1*d$eps_fam + gamma_2*rnorm(N)
    d$Y2 <- a + b*d$g2 + gamma_1*d$eps_fam + gamma_2*rnorm(N)
  } else {
    d$p0 <- sigmoid(a + b*d$g0 + gamma_1*d$eps_fam + gamma_2*rnorm(N))
    d$p1 <- sigmoid(a + b*d$g1 + gamma_1*d$eps_fam + gamma_2*rnorm(N))
    d$p2 <- sigmoid(a + b*d$g2 + gamma_1*d$eps_fam + gamma_2*rnorm(N))

    d$Y0 <- rbinom(N, 1, d$p0)
    d$Y1 <- rbinom(N, 1, d$p1)
    d$Y2 <- rbinom(N, 1, d$p2)
  }
  return(d)
}

r.sq <- function(p, b, gamma_1, gamma_2){
  top <- b^2 * 2*p*(1-p)
  bot <- top + gamma_1^2 + gamma_2^2
  r.sq <- top/bot
  return(r.sq)
}

sim.reg <- function(b, N, p, gamma_1, gamma_2){
  d <- sim.data(N=N, p_allele=p, continuous=TRUE,
                beta=b, gamma_1=gamma_1, gamma_2=gamma_2)
  d$g_pd <- 0.5 * d$g0
  d$delta_g1 <- d$g0 - d$g1
  d$delta_Y1 <- d$Y0 - d$Y1

  m_pd <- lm(delta_Y1 ~ g_pd, data=d)
  m_fe <- lm(delta_Y1 ~ delta_g1, data=d)
  m_ols <- lm(Y0 ~ g0, data=d)

  out <- c(summary(m_pd)$coefficients[2,4],
           coef(m_pd)[2],
           coef(m_fe)[2],
           coef(m_ols)[2])

  return(out)
}

# Simulation control parameters

p_allele <- 0.3
M <- 1000
b_len <- 50
bs_equal <- seq(from=0.02, to=0.2, length.out=b_len)
bs_unequal <- seq(from=0.03, to=0.3, length.out=b_len)
Ns <- c(5e2, 1e3, 5e3, 1e4, 5e4, 1e5)
N_len <- length(Ns)

gamma_1 <- c(1, 1, sqrt(3))
gamma_2 <- c(1, sqrt(3), 1)
case <- c('Equal', 'Individual', "Family")
N_cases <- length(case)

control <- data.frame(
  N = rep(Ns, each=M*b_len, times=N_cases),
  b = c(rep(bs_equal, each=M, times=N_len),
        rep(bs_unequal, each=M, times=(N_cases-1)*N_len)),
  p = rep(p_allele, times=b_len*N_len*M*N_cases),
  gamma_1 = rep(gamma_1, each=b_len*N_len*M),
  gamma_2 = rep(gamma_2, each=b_len*N_len*M),
  case = rep(case, each=b_len*N_len*M)
)


# Run simulation in parallel

cl <- makeCluster(20)
clusterSetRNGStream(cl, 1337)
clusterExport(cl, c('sim.reg', 'sim.data'))
sim_out <- clusterMap(cl, sim.reg,
                      b=control$b, N=control$N, p=control$p,
                      gamma_1=control$gamma_1, gamma_2=control$gamma_2,
                      SIMPLIFY=TRUE)
stopCluster(cl)

# Construct dataframe of results

results <- control
results$r_sq <- r.sq(results$p, results$b, results$gamma_1, results$gamma_2)
results$b_pd <- sim_out[2,]
results$b_fe <- sim_out[3,]
results$b_ols <- sim_out[4,]
results$p_val <- sim_out[1,]

saveRDS(results, '~/projects/phen-diff/sim_data.rds')
