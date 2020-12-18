library(tidyverse)

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }

  if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


se.ratio <- function(phi){
  r <- sqrt((1-phi)/(4-phi))
  return(r)
}

var.ratio <- function(phi){
  r <- (1-phi)/(4-phi)
  return(r)
}

pd.samp.se <- function(N_FE, K, phi){
  phi_v <- rep(phi, length(N_FE))
  N_PD <- sqrt(4-phi_v)*(N_FE-K-1) / (sqrt(1-phi_v)) + K + 1
  return(N_PD)
}

pd.samp.var <- function(N_FE, K, phi){
  phi_v <- rep(phi, length(N_FE))
  N_PD <- (4-phi_v)*(N_FE-K-1) / ((1-phi_v)) + K + 1
  return(N_PD)
}

phi <- seq(0,1, 0.01)

d <- tibble(
  phi = phi,
  ratio = se.ratio(phi))

N_fe <- seq(1, 50001, 100)

d2 <- tibble(
  N_FE = rep(N_fe, 4),
  Phi = as.factor(c(rep(0.01, length(N_fe)),
          rep(0.1, length(N_fe)),
          rep(0.3, length(N_fe)),
          rep(0.5, length(N_fe)))),
  N_PD = c(pd.samp.se(N_fe, 0, 0.01),
           pd.samp.se(N_fe, 0, 0.1),
           pd.samp.se(N_fe, 0, 0.3),
           pd.samp.se(N_fe, 0, 0.5))
)

p1 <- ggplot(d, aes(x = phi, y = ratio)) +
  geom_line(color='red') +
  labs(x = 'Phi',
       y = 'FE/PD Estimator SE Ratio',
       title = 'Estimator Standard Error Ratio by Phi') +
  coord_cartesian(xlim = c(0,1), ylim = c(0,0.5)) +
  theme_bw()

p2 <- ggplot(d2, aes(x = N_PD, y=N_FE, color=Phi)) +
  geom_line() +
  geom_abline(aes(intercept=0, slope=0.5), lty=2, alpha=0.4) +
  coord_cartesian(xlim = c(0,50000), ylim = c(0, 25000)) +
  labs(x = 'N PD Sibling Pairs',
       y = 'N FE Sibling Pairs',
       title = 'FE Sample Size Required to Match PD Precision') +
  theme_bw()

multiplot(p1, p2, cols=2)
