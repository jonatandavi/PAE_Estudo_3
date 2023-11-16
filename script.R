#--------------------------------------------------------

suppressPackageStartupMessages(library(smoof))

fn <- function(X){
  if(!is.matrix(X)) X <- matrix(X, nrow = 1) # <- if a single vector is passed as X
  
  Y <- apply(X, MARGIN = 1,
             FUN = smoof::makeRosenbrockFunction(dimensions = dim))
  
  return(Y)
  
}

## Config 1
recpars1 <- list(name = "recombination_exp", cr = 0.6)
mutpars1 <- list(name = "mutation_best", f = 2)

## Config 2
recpars2 <- list(name = "recombination_geo", alpha = 0.6)
mutpars2 <- list(name = "mutation_rand", f = 1.2)

set.seed(2)

dimensions <- sample(2:150, 64, replace=FALSE)
results1 <- matrix(0, 30, 64)

i <- 1
for(d in 1:64) {
  for(t in 1:30){
    dim <- dimensions[d]# d
    selpars <- list(name = "selection_standard")
    stopcrit <- list(names = "stop_maxeval", maxevals = 5000 * dim, maxiter = 100 * dim)
    probpars <- list(name = "fn", xmin = rep(-5, dim), xmax = rep(10, dim))
    popsize = 5 * dim
    
    suppressPackageStartupMessages(library(ExpDE))
    
    # Run algorithm on problem:
    out <- ExpDE(mutpars = mutpars1,
                 recpars = recpars1,
                 popsize = popsize,
                 selpars = selpars,
                 stopcrit = stopcrit,
                 probpars = probpars,
                 showpars = list(show.iters = "dots", showevery = 20))
    
    # Extract observation:
    out$Fbest
    results1[t, d] = out$Fbest
  }
  i = i + 1
}

results2 <- matrix(0, 30, 64)
i <- 1
print(i)
for(d in 1:64) {
  for(t in 1:30){
    dim <- dimensions[d]# d
    
    selpars <- list(name = "selection_standard")
    stopcrit <- list(names = "stop_maxeval", maxevals = 5000 * dim, maxiter = 100 * dim)
    probpars <- list(name = "fn", xmin = rep(-5, dim), xmax = rep(10, dim))
    popsize = 5 * dim
    
    suppressPackageStartupMessages(library(ExpDE))
    
    # Run algorithm on problem:
    out <- ExpDE(mutpars = mutpars2,
                 recpars = recpars2,
                 popsize = popsize,
                 selpars = selpars,
                 stopcrit = stopcrit,
                 probpars = probpars,
                 showpars = list(show.iters = "dots", showevery = 20))
    
    # Extract observation:
    out$Fbest
    
    results2[t, i] = out$Fbest
  }
  i = i + 1
}