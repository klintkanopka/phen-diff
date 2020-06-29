sim.phen <- function(N=1e5, alpha=0, beta=1, rho=0.5, 
                     mu_pgs=0, sd_pgs=1, sd_ind=1, sd_fam=1){
  
  sim.pgs <- function(N=1e5, rho=0.5, mu=0, sd=1){
    require(MASS)
    var_mat <- matrix(c(sd^2, rho*sd^2, rho*sd^2, sd^2), ncol=2)
    pgs <- mvrnorm(n=N, mu=c(mu,mu), Sigma=var_mat, empirical=TRUE)
    out <- data.frame(g0 = c(pgs[,1], pgs[,2]),
                      g1 = c(pgs[,2], pgs[,1]),
                      family = c(1:N, 1:N),
                      primary = c(rep(1, N), rep(0, N)))
    return(out)
  }
  
  eps_0_tmp <- rnorm(N, sd=sd_ind)
  eps_1_tmp <- rnorm(N, sd=sd_ind)
  eps_fam_tmp <- rnorm(N, sd=sd_fam)
  
  eps_0 <- c(eps_0_tmp, eps_1_tmp)
  eps_1 <- c(eps_1_tmp, eps_0_tmp)
  eps_fam <- rep(eps_fam_tmp, 2)
  a <- rep(alpha, 2*N)
  b <- rep(beta, 2*N)
  
  d <- sim.pgs(N=N, rho=rho, mu=mu_pgs, sd=sd_pgs)
  
  d$Y0 <- a + b*d$g0 + eps_0 + eps_fam
  d$Y1 <- a + b*d$g1 + eps_1 + eps_fam
  
  return(d)
}