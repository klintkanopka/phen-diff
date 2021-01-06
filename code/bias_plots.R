library(tidyverse)

results <- readRDS('code/sim_data/sim_allele_cont.rds')

results %>%
  filter(N==500) %>%
  select(N, b, case, b_pd, b_fe, b_ols) %>%
  pivot_longer(starts_with('b_'), names_to='model', values_to='estimate') %>%
  ggplot(aes(x=b, y=estimate)) +
  geom_point(color='grey', alpha=0.05) +
  geom_smooth(aes(color=model), method=lm, formula=y~x) +
  geom_abline(aes(slope=1, intercept=0), lty=2, color='black') +
  facet_grid(model~case) +
  theme_bw()

ggsave('figs/allele_bias_500.png', width=8, height=6)
