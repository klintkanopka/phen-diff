# simulates a continuous phenotype that is related to a genotype expressed
# as an allele count for three siblings suitable for OLS, FE, or phenotype
# differences estimation
#
# inputs:
#   N - number of sibling triplets
#   rho_g - genetic correlation in parent's generation (approximate)
#   p_allele - allele prevalence in parent's generation
#   alpha - mean phenotype for genotype 0
#   beta - magnitude of direct genetic effect
#   gamma_1 - standard deviation of family environment effect distribution
#   gamma_2 - standard deviation of individual environment effect distribution
#
# outputs:
#   dataframe with:
#     g_m/g_p - maternal and paternal genotypes
#     g0/g1/g2 - sibling genotypes
#     eps_fam - family environment
#     Y0/Y1/Y2 - sibling phenotypes


sim.allele.cont <- function(N=1e5, rho_g=0, p_allele=0.3,
                            alpha=0, beta=0.04,
                            gamma_1=sqrt(0.35), gamma_2=sqrt(0.649)){

  sim.genes <- function(N=N, rho=rho, rho_g=rho_g, p_allele=p_allele){

    g_m <- rbinom(n=N, size=2, prob=p_allele)
    g_p <- rbinom(n=N, size=2, prob=(1-rho_g)*p_allele + rho_g*(g_m/2))

    out <- data.frame(g_m = g_m,
                      g_p = g_p,
                      g0 = rbinom(N, 1, g_m/2) + rbinom(N, 1, g_p/2),
                      g1 = rbinom(N, 1, g_m/2) + rbinom(N, 1, g_p/2),
                      g2 = rbinom(N, 1, g_m/2) + rbinom(N, 1, g_p/2),
                      eps_fam = rnorm(N))
    return(out)
  }

  d <- sim.genes(N=N, rho=rho, rho_g=rho_g, p_allele=p_allele)

  d$Y0 <- alpha + beta*d$g0 + gamma_1*d$eps_fam + gamma_2*rnorm(N)
  d$Y1 <- alpha + beta*d$g1 + gamma_1*d$eps_fam + gamma_2*rnorm(N)
  d$Y2 <- alpha + beta*d$g2 + gamma_1*d$eps_fam + gamma_2*rnorm(N)

  return(d)
}
