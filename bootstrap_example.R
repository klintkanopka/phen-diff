library(MASS)
library(tidyverse)

source('sim_utils.R')

# compute bootstrap SEs - these are the default parameters
out <- bootstrap.sim(M=1e3, N=1e2, # bootstrap iterations and observations
                     alpha=0, beta=1, gamma=1, # params of causal model
                     rho=0.5, rho_g=0.7, # sibling and genetic nurture corrs
                     mu_dir=0, sd_dir=1, mu_nur=0, sd_nur=NULL, # effect dist
                     sd_ind=1, sd_fam=1, # error params
                     het = 'none' # heteroskedasticity form
)

# the out object contains the betas, a summary of the ses, and a plot
out$betas
out$summary
out$plot

# we can try more iterations, obeservations, and heteroskedasticity 
out <- bootstrap.sim(M=1e4, N=1e3, het='full')
# you can also add stuff to the plot
out$plot + 
  labs(x='Beta Hat', y='Density', title='PD is less efficient') +
  scale_fill_manual(values=c('FE' = 'cornflowerblue',
                             'PD' = 'darkgreen',
                             'PD_full' = 'darkorchid')) +
  scale_color_manual(values=c('FE' = 'cornflowerblue',
                              'PD' = 'darkgreen',
                              'PD_full' = 'darkorchid'))
