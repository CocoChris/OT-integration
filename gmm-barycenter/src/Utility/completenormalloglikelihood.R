complete.normal.loglikelihood <- function(x, z, pars){
  n <- length(x)
  logl =.loopexp_scalar_completeL(log(pars[ , 3]), pars[ , 1],
                                  sqrt(pars[ , 2]), z, x,n)
  # g <- dim(pars)[1]
  # 
  # logl <- rep(0, n)
  # logpi <- log(pars[ , 3])
  # mean <- pars[ , 1]
  # sigma <- sqrt(pars[ , 2])
  # logl <- logpi[z] + dnorm(x, mean = mean[z], sd = sigma[z], log = T)
  return(sum(logl))}