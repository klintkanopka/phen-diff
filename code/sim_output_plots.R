library(tidyverse)

results <- readRDS('sim_data.rds')

results %>%
  filter(N>5000) %>%
  ggplot(aes(x=b, y=b_pd)) +
  geom_point(color='grey', alpha=0.05) +
  geom_smooth(color='dodgerblue1', method=lm, formula=y~x) +
  geom_abline(aes(slope=1, intercept=0), lty=2, color='black') +
  facet_grid(case~.) +
  theme_bw()

results %>%
  group_by(N, b, case) %>%
  summarize(b_hat = mean(b_pd),
            .groups='drop') %>%
  mutate(error = (b_hat-b)/b) %>%
  ggplot(aes(x=b, y=error)) +
  geom_point(color='grey', alpha=0.75) +
  geom_smooth(method=lm, formula=y~x) +
  facet_grid(.~case) +
  theme_bw()
