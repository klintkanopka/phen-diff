library(tidyverse)

results <- readRDS('code/sim_data/sim_allele_bias.rds')

results %>%
  select(N, b, case, b_pd, b_fe, b_ols) %>%
  pivot_longer(starts_with('b_'), names_to='model', values_to='estimate') %>%
  mutate(case = factor(case, levels=c('None', 'Low', 'High')),
         mod = case_when(model == 'b_fe' ~ 'Fixed Effects',
                         model == 'b_ols' ~ 'OLS',
                         model == 'b_pd' ~ 'Phenotype Differences'),
         mod = factor(mod, levels=c('OLS', 'Fixed Effects', 'Phenotype Differences'))) %>%
  filter(case == 'High') %>%
  ggplot(aes(x=b, y=estimate, color=mod)) +
  geom_point(alpha=0.015) +
  geom_smooth(method=lm, formula=y~x, se=F) +
  geom_abline(aes(slope=1, intercept=0), lty=2, color='black', alpha=0.3) +
  facet_grid(.~mod) +
  labs(x='Simulated Beta', y='Estimate',
       title='Bias in Beta estimates by level of gene-environment correlation and estimator',
       color='Estimator') +
  theme_bw() +
  theme(legend.position="bottom")

ggsave('figs/high_allele_bias.png', width=8, height=6)
