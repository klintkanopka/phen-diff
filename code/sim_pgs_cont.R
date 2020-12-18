# simulates a continuous phenotype that is related to a genotype expressed
# as a polygenic score for three siblings suitable for OLS, FE, or phenotype
# differences estimation
#
# inputs:
#   N - number of sibling triplets
#   rho - sibling genetic correlation
#   rho_g - gene-environment correlation
#   alpha - mean phenotype for genotype 0
#   beta - magnitude of direct genetic effect
#   gamma_1 - standard deviation of family environment effect distribution
#   gamma_2 - standard deviation of individual environment effect distribution
#   gamma_3 - standard deviation of genetic nurture effect
#
# outputs:
#   dataframe with:
#     g0/g1/g2 - sibling genotypes
#     g_nur - magnitude of family genetic nurture effect
#     eps_fam - family environment
#     Y0/Y1/Y2 - sibling phenotypes


sim.pgs.cont <- function(N=1e5, rho=0.5, rho_g=0.7,
                         alpha=-0.85, beta=1,
                         gamma_1=0.5, gamma_2=0.25, gamma_3=NULL){

  sim.pgs <- function(N=1e5, rho=0.5, rho_g=0.7, gamma_3=NULL){
    require(MASS)
    if (is.null(gamma_3)){gamma_3 <- 0.7}
    var_mat <- matrix(c(1, rho, rho, rho_g*gamma_3,
                        rho, 1, rho, rho_g*gamma_3,
                        rho, rho, 1, rho_g*gamma_3,
                        rho_g*gamma_3, rho_g*gamma_3, rho_g*gamma_3, gamma_3^2),
                      ncol=4)
    pgs <- mvrnorm(n=N, mu=c(0, 0, 0, 0), Sigma=var_mat, empirical=TRUE)
    out <- data.frame(g0 = pgs[,1],
                      g1 = pgs[,2],
                      g2 = pgs[,3],
                      g_nur = pgs[,4],
                      eps_fam = rnorm(N))
    return(out)
  }

  eps_fam <- rnorm(N)

  d <- sim.pgs(N=N, rho=rho, rho_g=rho_g, gamma_3=gamma_3)

  d$Y0 <- alpha + beta*d$g0 + gamma_1*d$eps_fam + gamma_2*rnorm(N) + d$g_nur
  d$Y1 <- alpha + beta*d$g1 + gamma_1*d$eps_fam + gamma_2*rnorm(N) + d$g_nur
  d$Y2 <- alpha + beta*d$g2 + gamma_1*d$eps_fam + gamma_2*rnorm(N) + d$g_nur

  return(d)
}
