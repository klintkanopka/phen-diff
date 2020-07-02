sim.phen <- function(N=1e5, alpha=0, beta=1, rho=0.5, rho_g=0.7, gamma=1, 
                     mu_dir=0, sd_dir=1, mu_nur=0, sd_nur=NULL, 
                     sd_ind=1, sd_fam=1,
                     het = 'none'){
  
  sim.pgs <- function(N=1e5, rho=0.5, rho_g=0.7, 
                      mu_dir=0, sd_dir=1, 
                      mu_nur=0, sd_nur=NULL){
    require(MASS)
    if (is.null(sd_nur)){sd_nur <- 0.7*sd_dir}
    var_mat <- matrix(c(sd_dir^2, rho*sd_dir^2, rho_g*sd_dir*sd_nur, 
                        rho*sd_dir^2, sd_dir^2, rho_g*sd_dir*sd_nur,
                        rho_g*sd_dir*sd_nur, rho_g*sd_dir*sd_nur, sd_nur^2), 
                      ncol=3)
    pgs <- mvrnorm(n=N, mu=c(mu_dir,mu_dir,mu_nur), 
                   Sigma=var_mat, empirical=TRUE)
    out <- data.frame(g_d0 = c(pgs[,1], pgs[,2]),
                      g_d1 = c(pgs[,2], pgs[,1]),
                      g_nur = rep(pgs[,3], 2),
                      family = rep(1:N, 2),
                      primary = c(rep(1, N), rep(0, N)))
    return(out)
  }

  if (het %in% c('full', 'individual')){
    ind <- 1
  } else {
    ind <- 0
  }
  
  if (het %in% c('full', 'family')){
    fam <- 1
  } else {
    fam <- 0
  }

  eps_0_tmp <- rnorm(N, sd=sd_ind)
  eps_1_tmp <- rnorm(N, sd=sd_ind)
  eps_fam_tmp <- rnorm(N, sd=sd_fam)
  
  eps_0 <- c(eps_0_tmp, eps_1_tmp)
  eps_1 <- c(eps_1_tmp, eps_0_tmp)
  eps_fam <- rep(eps_fam_tmp, 2)
  a <- rep(alpha, 2*N)
  b <- rep(beta, 2*N)
  
  d <- sim.pgs(N=N, rho=rho, rho_g=rho_g, 
               mu_dir=mu_dir, sd_dir=sd_dir, 
               mu_nur=mu_nur, sd_nur=sd_nur)
  
  d$Y0 <- a + b*d$g_d0 + gamma*d$g_nur + 
    eps_0*(1+ind*d$g_d0) + eps_fam*(1+fam*d$g_d0)
  d$Y1 <- a + b*d$g_d1 + gamma*d$g_nur + 
    eps_1*(1+ind*d$g_d1) + eps_fam*(1+fam*d$g_d1)
  
  return(d)
}


bootstrap.sim <- function(M=1e3, N=1e2, 
                          alpha=0, beta=1, gamma=1,
                          rho=0.5, rho_g=0.7,
                          mu_dir=0, sd_dir=1, mu_nur=0, sd_nur=NULL, 
                          sd_ind=1, sd_fam=1,
                          het = 'none'
                          ){
  require(tidyverse)
  
  d <- sim.phen(N=N)
  
  d$delta_g <- d$g_d1 - d$g_d0
  d$delta_Y <- d$Y1 - d$Y0
  
  d_lim <- d[d$primary == 1,]
  
  betas_FE <- c()
  betas_PD <- c()
  betas_PD_full <- c()
  
  for (i in 1:M){
    
    samp <- sample(1:N, N, replace=TRUE)
    samp_full <- c(samp, samp+N)
    
    boot_FE <- lm(delta_Y ~ delta_g, data=d[samp,])
    boot_PD <- lm(delta_Y ~ g_d1, data=d[samp,])
    boot_PD_full <- lm(delta_Y ~ g_d1, data=d[samp_full,])
    
    betas_FE[i] <- coef(boot_FE)[2]
    betas_PD[i] <- coef(boot_PD)[2]
    betas_PD_full[i] <- coef(boot_PD_full)[2]
  }
  
  betas <- data.frame(FE = betas_FE,
                      PD = betas_PD,
                      PD_full = betas_PD_full)
  
  tmp <- betas %>% 
    mutate(PD = 2*PD,
           PD_full = 2*PD_full) %>% 
    pivot_longer(cols=everything(), names_to = 'model', values_to='beta') 
  
  summary <- tmp %>%
    group_by(model) %>% 
    summarize(beta_hat = mean(beta), 
              se = sd(beta),
              .groups='drop')
  
  plot <- ggplot(tmp, aes(x = beta, color = model, fill = model)) +
    geom_density(alpha = 0.6) +
    theme_bw()
  print(summary)
  out <- list(betas=betas, summary=summary, plot=plot)
  return(out)
}