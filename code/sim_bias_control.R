library(tidyverse)
library(parallel)

source('code/sim_allele_cont.R')
source('code/r_sq.R')

sim.reg <- function(b, N, p, gamma_1, gamma_2, rho_g){
  d <- sim.allele.cont(N=N,p_allele=p,beta=b,gamma_1=gamma_1,gamma_2=gamma_2,rho_g=rho_g)
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
rho_gs <- c(0, 0.2, 0.5)
M <- 1000
b_len <- 50
bs <- seq(from=0.02, to=0.2, length.out=b_len)
Ns <- c(1e3)
N_len <- length(Ns)

gamma_1 <- c(1)
gamma_2 <- c(1)
case <- c('None', 'Low', "High")
N_cases <- length(case)

control <- tibble(
  rho_g = rep(rho_gs, each=M*b_len),
  N = 500,
  b = c(rep(bs, each=M, times=N_cases)),
  p = p_allele,
  gamma_1 = 1,
  gamma_2 = 1,
  case = rep(case, each=b_len*M)
)


# Run simulation in parallel

cl <- makeCluster(20)
clusterSetRNGStream(cl, 1337)
clusterExport(cl, c('sim.reg', 'sim.allele.cont'))
sim_out <- clusterMap(cl, sim.reg, rho_g=control$rho_g,
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

saveRDS(results, 'code/sim_data/sim_allele_bias.rds')
