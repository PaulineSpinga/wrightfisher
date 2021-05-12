popSimu <- function(N, nbA = N) # par défaut nbA = N
{
  N_tot <- 2*N
  if(N_tot == 0){stop('N_tot equals 0')}
  if(nbA > N_tot){stop('nbA greater than N_tot')}
  if(nbA < 0){stop('nbA is negative')}
  if (!(is.numeric(nbA) && (is.numeric(N_tot)))){stop('one argument is not a numeric')}


  p <- (nbA/(N_tot))
  vectA <- nbA
  h <- (2*N * p * (1-p))/(2*N-1)
  while((nbA > 0) && (nbA < (2*N)))
  {
    #new Allele repartition
    pop <- rbinom(n=2*N, size=1, prob=p) #A = 1, a = 0
    #update nbA and p
    nbA <- sum(pop)
    p <- nbA/(2*N) #  p <- sum(pop)/(2*N)
    #save result in vectA and h
    vectA <- c(vectA, nbA)
    h <- c(h, (2*N * p * (1-p))/(2*N-1))
  }
  return (list(nbA = vectA, fixationTime = length(vectA), coeffH = h))
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

createData <- function(tirages)
{
  nb_row <- maxA(tirages)
  nb_tir <- length(tirages)

  df <- data.frame(y=rep(1,nb_row))
  for(i in 1:nb_tir)
  {
    temp_df <- data.frame(y=rep(0,nb_row))
    row_i <- tirages[[i]]$fixationTime
    temp_df[1:row_i,] <- tirages[[i]]$nbA
    df <- cbind(df,temp_df)
    df[,i] <- temp_df
  }
  return(df)
}

info_T <- function(vect_temps){
  esp <- mean(vect_temps)
  variance <- var(vect_temps)

  return(list(esperance=esp, variance=variance))

}

mclapply.simu <- function(Ne,Na, nbRep){
  listout <- vector("list", nbRep)
  for(i in 1:nbRep){
    listout[[i]] <- popSimu(Ne,Na)
  }
  return(listout)
}


