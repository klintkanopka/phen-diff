library(tidyverse)
library(parallel)

source('code/sim_allele_cont.R')
source('code/r_sq.R')

sim.reg.bias.error <- function(b, N, p, gamma_1, gamma_2, rho_g, bias, sd_error){
  d <- sim.allele.cont(N=N,p_allele=p,beta=b,gamma_1=gamma_1,gamma_2=gamma_2,rho_g=rho_g)
  d$g_pd <- 0.5 * d$g0
  error <- rnorm(N, mean=0, sd=sd_error)
  d$delta_Y1 <- d$Y0 - d$Y1
  d$delta_Y1_e <- d$Y0 - (d$Y1 + error)
  d$delta_Y1_b <- d$Y0 - (d$Y1 + bias)
  d$delta_Y1_be <- d$Y0 - (d$Y1 + error + bias)

  m_pd <- lm(delta_Y1 ~ g_pd, data=d)
  m_pd_e <- lm(delta_Y1_e ~ g_pd, data=d)
  m_pd_b <- lm(delta_Y1_b ~ g_pd, data=d)
  m_pd_be <- lm(delta_Y1_be ~ g_pd, data=d)

  out <- c(coef(m_pd)[2],
           coef(m_pd_e)[2],
           coef(m_pd_b)[2],
           coef(m_pd_be)[2])

  return(out)
}

# Simulation control parameters

M <- 5000
bias = seq(0.1, 2, by=0.1)
sd_error = seq(0.1, 2, by=0.1)

control <- tibble(
  rho_g = 0.3,
  N = 500,
  b = 0.05,
  p = 0.3,
  gamma_1 = 1,
  gamma_2 = 1,
  bias = rep(bias, each=M*length(bias)),
  sd_error = rep(sd_error, times=M*length(sd_error))
)


# Run simulation in parallel

cl <- makeCluster(20)
clusterSetRNGStream(cl, 1337)
clusterExport(cl, c('sim.reg.bias.error', 'sim.allele.cont'))
sim_out <- clusterMap(cl, sim.reg.bias.error, rho_g=control$rho_g,
                      b=control$b, N=control$N, p=control$p,
                      gamma_1=control$gamma_1, gamma_2=control$gamma_2,
                      bias=control$bias, sd_error=control$sd_error,
                      SIMPLIFY=TRUE)
stopCluster(cl)

# Construct dataframe of results

results <- control
results$r_sq <- r.sq(results$p, results$b, results$gamma_1, results$gamma_2)
results$b_pd <- sim_out[1,]
results$b_pd_e <- sim_out[2,]
results$b_pd_b <- sim_out[3,]
results$b_pd_be <- sim_out[4,]


saveRDS(results, 'code/sim_data/sim_allele_bias_error.rds')


ggplot(results, aes(x = jitter(bias), y = b_pd_b)) +
  geom_point(alpha = 0.005, color='dodgerblue1') +
  geom_smooth(method=lm, formula=y~x, color='firebrick1', se=F) +
  geom_hline(aes(yintercept=0.05), lty=2, color='black') +
  theme_bw()

ggsave('figs/unbiased_bias.png', width=8, height=6)


ggplot(results, aes(x = jitter(sd_error), y = b_pd_e)) +
  geom_point(alpha = 0.005, color='dodgerblue1') +
  geom_smooth(method=lm, formula=y~x, color='firebrick1', se=F) +
  geom_hline(aes(yintercept=0.05), lty=2, color='black') +
  theme_bw()

ggsave('figs/unbiased_error.png', width=8, height=6)


results %>%
  group_by(bias, sd_error) %>%
  summarize(mean = mean(b_pd_be),
            sd = sd(b_pd_be)) %>%
  ggplot(aes(x=bias, y=sd_error, fill=mean)) +
  geom_tile() +
  scale_fill_gradient2(low='dodgerblue1', mid='firebrick1', high='goldenrod1',
                       midpoint=0.05) +
  labs(x = 'Bias', y = 'Error SD', fill='Mean Beta',
       title='Unbiased under both conditions (N=5000 per cell)',
       subtitle='Target value is red (0.05)') +
  theme_bw()

ggsave('figs/unbiased_bias_error.png', width=8, height=6)


results %>%
  group_by(bias, sd_error) %>%
  summarize(mean = mean(b_pd_be),
            sd = sd(b_pd_be)) %>%
  ggplot(aes(x=bias, y=sd_error, fill=sd)) +
  geom_tile() +
  scale_fill_gradient2(low='dodgerblue1', mid='firebrick1', high='goldenrod1',
                       midpoint=0.275) +
  labs(x = 'Bias', y = 'Error SD', fill='Beta SE',
       title='SE on Beta estimate increases with error SD, not bias (N=5000)') +
  theme_bw()

ggsave('figs/standard_error_bias_error.png', width=8, height=6)

