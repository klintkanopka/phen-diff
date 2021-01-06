r.sq <- function(p, b, gamma_1, gamma_2){
  top <- b^2 * 2*p*(1-p)
  bot <- top + gamma_1^2 + gamma_2^2
  r.sq <- top/bot
  return(r.sq)
}
