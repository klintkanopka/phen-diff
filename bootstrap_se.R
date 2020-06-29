library(tidyverse)

setwd('~/projects/phen-diff')
source('sim_data.R')

N <- 1e4
M <- 1e4

d <- sim.phen(N=N)

d$delta_g <- d$g1 - d$g0
d$delta_Y <- d$Y1 - d$Y0

d_lim <- d[d$primary == 1, ]

m_FE <- lm(delta_Y ~ delta_g, data=d_lim)
m_PD <- lm(delta_Y ~ g1, data=d_lim)
m_PD_full <- lm(delta_Y ~ g1, data=d)

betas_FE <- c()
betas_PD <- c()
betas_PD_full <- c()
 
for (i in 1:M){
  
  samp <- sample(1:N, N, replace=TRUE)
  samp_full <- c(samp, samp+N)
  
  boot_FE <- lm(delta_Y ~ delta_g, data=d[samp,])
  boot_PD <- lm(delta_Y ~ g1, data=d[samp,])
  boot_PD_full <- lm(delta_Y ~ g1, data=d[samp_full,])
  
  betas_FE[i] <- coef(boot_FE)[2]
  betas_PD[i] <- coef(boot_PD)[2]
  betas_PD_full[i] <- coef(boot_PD_full)[2]
}

betas <- data.frame(FE = betas_FE,
                    PD = betas_PD,
                    PD_full = betas_PD_full)

summary(betas)

betas %>% 
  mutate(PD = 2*PD,
         PD_full = 2*PD_full) %>% 
  pivot_longer(cols=everything(), names_to = 'model', values_to='beta') %>%
  group_by(model) %>% 
  summarize(beta_hat = mean(beta), 
                se = sd(beta))
