relabelmcmc <- function(set, m, burn,mcmc.pars, ls){
  ll <- length(set)
  mcmcsample <- 1:(m-burn)
  output <- list()

  for(j in 1:ll){
    assign(paste0("mcmc.pars.", set[j], sep = ""), mcmc.pars[mcmcsample,,])
    dv <- get(paste0("mcmc.pars.", set[j], sep = ""))
    eval(parse(text=paste0("rm(", "mcmc.pars.",set[j],")", sep = "")))
    # remove the original object from memory
    lset <- which(names(ls$permutations) == set[j])
    dset <- ls$permutations[[lset]]
    for(it in 1:(m-burn)){
      dv[it,,] <- dv[it,dset[mcmcsample[it],],]
    }
    assign(paste("mcmc.pars.", set[j], sep = ""), dv)
    output[[paste0("mcmc.pars.", set[j], sep = "")]] <- dv
  }
  return(output)
}