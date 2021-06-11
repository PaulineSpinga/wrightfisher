popSimu <- function(N, nbA = N)
{
  N_tot <- 2*N
  if(N_tot == 0){stop('N_tot equals 0')}
  if(nbA > N_tot){stop('nbA greater than N_tot')}
  if(nbA < 0){stop('nbA is negative')}
  if (!(is.numeric(nbA) && (is.numeric(N_tot)))){stop('one argument is not a numeric')}

  proba_init <- (nbA/(N_tot))
  p <- proba_init
  vectA <- nbA
  h <- (2 * p * (1-p))
  while((nbA > 0) && (nbA < (2*N)))
  {
    pop <- rbinom(n=2*N, size=1, prob=p) #A = 1, a = 0
    nbA <- sum(pop)
    p <- nbA/(2*N)
    vectA <- c(vectA, nbA)
    h <- c(h, (2 * p * (1-p)))
  }

  temps <- length(vectA)-1

  return (list(nbA = vectA, fixationTime = temps, coeffH = h, proba_init = proba_init))
}


popSimu.select <- function(N, nbA = N,s) # par défaut nbA = N
{
  N_tot <- 2*N
  if(N_tot == 0){stop('N_tot equals 0')}
  if(nbA > N_tot){stop('nbA greater than N_tot')}
  if(nbA < 0){stop('nbA is negative')}
  if (!(is.numeric(nbA) && (is.numeric(N_tot)))){stop('one argument is not a numeric')}

  proba_init <- (nbA/(N_tot))
  p <- proba_init
  vectA <- nbA
  h <- (2 * p * (1-p))
  while((nbA > 0) && (nbA < (2*N)))
  {
    pop <- rbinom(n=2*N, size=1, prob=p) #A = 1, a = 0
    nbA <- sum(pop)
    p <- (1+s)*nbA/((1+s)*nbA + N_tot - nbA)
    vectA <- c(vectA, nbA)
    h <- c(h, (2 * p * (1-p)))
  }

  temps <- length(vectA) - 1

  return (list(nbA = vectA, fixationTime = temps, coeffH = h, proba_init = proba_init))
}



maxA <- function(tirages)
{
  max_A <- NULL
  for (i in 1:length(tirages))
  {
    max_A <- c(max_A,tirages[[i]]$fixationTime)
  }
  return(max(max_A))
}

vectFix <- function(tirages)
{
  vect <- NULL
  for (i in 1:length(tirages))
  {
    vect <- c(vect,tirages[[i]]$fixationTime)
  }
  return(vect)
}



vectPropa <- function(listSimu){
  vect <- NULL
  for (i in 1:length(listSimu)){
    vect <- c(vect, listSimu[[i]][[1]]$proba_init)
  }
  return(vect)
}



vectTempsMoy <- function(listSimu){
  vect <- NULL
  for (i in 1:length(listSimu)){
    vect <- c(vect,mean(vectFix(listSimu[[i]])))
  }
  return(vect)
}

simu.simu <- function (Ne, step, Nb_rep){
  if(step >= 2*Ne){stop('invalid step')}
  valA <- seq(0,2*Ne, by=step)
  outF <- vector("list", length(valA))
  for (i in 1:length(valA)){
    outF[[i]] <- mclapply(rep(Ne,Nb_rep),popSimu, nbA = valA[i],mc.cores = 1)
  }
  return(outF)
}

simu.simu.select <- function (Ne, step, Nb_rep,s){
  if(step >= 2*Ne){stop('invalid step')}
  valA <- seq(0,2*Ne, by=step)
  outF <- vector("list", length(valA))
  for (i in 1:length(valA)){
    outF[[i]] <- mclapply(rep(Ne,Nb_rep),popSimu.select, nbA = valA[i],s=s,mc.cores = 1)
  }
  return(outF)
}

