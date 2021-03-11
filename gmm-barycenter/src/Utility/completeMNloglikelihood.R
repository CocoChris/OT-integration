complete.md.normal.loglikelihood <- function(x, z, pars) {  
  # x: denotes the n data points
  # z: denotes an allocation vector (size=n)
  # pars: K\times 3 vector of means,variance, weights
  # pars[k,1]: corresponds to the mean of component k
  # pars[k,2]: corresponds to the variance of component k
  # pars[k,3]: corresponds to the weight of component k
  g <- dim(pars)[1]  ## number of components
  n <- dim(x)[1]
  dx <- dim(x)[2]
  J = dx + dx*(dx+1)/2 + 1
  logl <- rep(0, n)
  logpi <- log(pars[, J])
  my_mean <- pars[, 1:dx]
  logl <- logpi[z]
  for (k in 1:g) {
    ind <- which(z == k)
    if (length(ind) > 0) {
      cov.mat <- array(data = 0, dim = c(dx, dx))
      cov.mat[lower.tri(cov.mat, diag = T)] <- cov.mat[upper.tri(cov.mat, diag = T)] <-pars[k,(dx+1):(J-1)]
      logl[ind] <- logl[ind] + dmvnorm(x[ind, ], mean = my_mean[k, ], sigma = cov.mat, 
                                       log = TRUE)
    }
  }
  return(sum(logl))
}
