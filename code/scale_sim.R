library(tidyverse)
library(parallel)

source('sim_allele_cont.R')

r.sq <- function(p, b, gamma_1, gamma_2){
  top <- b^2 * 2*p*(1-p)
  bot <- top + gamma_1^2 + gamma_2^2
  r.sq <- top/bot
  return(r.sq)
}

sim.scale.reg <- function(b, N, p, gamma_1, gamma_2){
  d <- sim.allele.cont(N=N, p_allele=p, continuous=TRUE,
                       beta=b, gamma_1=gamma_1, gamma_2=gamma_2)
  d$delta_Y1 <- d$Y0 - d$Y1

  m_pd <- lm(delta_Y1 ~ g0, data=d)

  out <- c(summary(m_pd)$coefficients[2,4],
           coef(m_pd)[2])

  return(out)
}

# Simulation control parameters

p_allele <- seq(from=0.1, to=0.5, by=0.1)
M <- 1000
b_len <- 50
bs <- seq(from=0.02, to=0.2, length.out=b_len)
Ns <- c(1e3)
N_len <- length(Ns)

gamma_1 <- c(1)
gamma_2 <- c(1)
case <- c('Equal')
N_cases <- length(case)

control <- data.frame(
  N = rep(Ns, each=M*b_len, times=N_cases),
  b = rep(bs, each=M, times=N_len),
  p = rep(p_allele, times=b_len*N_len*M*N_cases),
  gamma_1 = rep(gamma_1, each=b_len*N_len*M),
  gamma_2 = rep(gamma_2, each=b_len*N_len*M),
  case = rep(case, each=b_len*N_len*M)
)


# Run simulation in parallel

cl <- makeCluster(20)
clusterSetRNGStream(cl, 1337)
clusterExport(cl, c('sim.scale.reg', 'sim.data'))
sim_out <- clusterMap(cl, sim.scale.reg,
                      b=control$b, N=control$N, p=control$p,
                      gamma_1=control$gamma_1, gamma_2=control$gamma_2,
                      SIMPLIFY=TRUE)
stopCluster(cl)

# Construct dataframe of results

results <- control
results$r_sq <- r.sq(results$p, results$b, results$gamma_1, results$gamma_2)
results$b_pd <- sim_out[2,]
results$p_val <- sim_out[1,]

saveRDS(results, '~/projects/phen-diff/sim_scale_data.rds')


# new sim: what's the correct scale factor?

results <- results %>%
  mutate(ratio = b_pd/b)


summary(lm(ratio~p, data=results))
summary(lm(ratio~b, data=results))


results %>%
  select(b, p, b_pd, ratio) %>%
  ggplot(aes(x=b, y=b_pd, color=as.factor(p))) +
  geom_point(color='grey', alpha=0.1) +
  geom_smooth(method=lm, formula=y~x) +
  geom_abline(aes(intercept=0, slope=1), color='black', lty=2) +
  theme_bw()


results %>%
  select(b, p, b_pd, ratio) %>%
  pivot_longer(c(b,p), names_to='measure', values_to='value') %>%
  ggplot(aes(x=value, y=ratio)) +
  geom_point(color='grey', alpha=0.1) +
  geom_smooth(method=lm, formula=y~x) +
  geom_hline(aes(yintercept=0.5), color='black', lty=2) +
  facet_grid(measure~.) +
  theme_bw()

results %>%
  select(b, p, ratio) %>%
  group_by(b, p) %>%
  summarize(mean_ratio = mean(ratio),
            .groups='drop') %>%
  pivot_longer(c(b,p), names_to='measure', values_to='value') %>%
  ggplot(aes(x=value, y=mean_ratio)) +
  geom_point(alpha=0.5) +
  geom_hline(aes(yintercept=0.5), color='black', lty=2) +
  facet_grid(measure~.) +
  theme_bw()

mean(results$ratio)
