#simulation du temps de fixation pour une population
#on tire une nouvelle génération tant p différent de 0 et 1
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


#pareil que popSimu mais avec ajout de l'effet de sélection
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


#fonction qui retourne le temps de fixation le plus grand (et donc le nombre max d'allèles A ajoutés)
#parmi toutes les simulations réalisées pour un jeu de paramètres (taille pop et nombre d'allèles A initiaux)
maxA <- function(tirages)
{
  max_A <- NULL
  for (i in 1:length(tirages))
  {
    max_A <- c(max_A,tirages[[i]]$fixationTime)
  }
  return(max(max_A))
}

##fonction qui retourne le vecteur contenant des temps de fixation pour chaque popSimu réalisée
vectFix <- function(tirages)
{
  vect <- NULL
  for (i in 1:length(tirages))
  {
    vect <- c(vect,tirages[[i]]$fixationTime)
  }
  return(vect)
}

##fonction qui renvoie de vecteur de proba initiales de chaque set de simulation
## 1 set de simulation : plusieurs popSimu (ou popSimu.select) réalisées pour plusieurs proba initiales différentes
vectPropa <- function(listSimu){
  vect <- NULL
  for (i in 1:length(listSimu)){
    vect <- c(vect, listSimu[[i]][[1]]$proba_init)
  }
  return(vect)
}

##fonction qui renvoie le temps moyen de chaque set de simu

vectTempsMoy <- function(listSimu){
  vect <- NULL
  for (i in 1:length(listSimu)){
    vect <- c(vect,mean(vectFix(listSimu[[i]])))
  }
  return(vect)
}

################################

#fonction non utilisée
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

############# ggplot #####################

#nombre d'allèles A en fonction du temps (pour popSimu et popSimu.select)
affichDataA <- function(tirages,Nb_rep){

  df <- NULL
  for(i in 1:Nb_rep)
  {
    len <- tirages[[i]]$fixationTime
    temp_df <- data.frame(x=0:len, y=tirages[[i]]$nbA, col=rep(i, each=len+1))
    df <- rbind(df,temp_df)
  }
  ggplot(df,aes(x=x,y=y,group=col,colour=factor(col))) + geom_line(size=1) + labs(title="Number of A over time",x="Time", y = "Number of A") +
    theme(axis.title=element_text(size=14,face="bold"),plot.title=element_text(size=18,face="bold",hjust=0.5),legend.position = "none")  # plot data

}

#coefficient d'hétérozygotie en fonction du temps (pour popSimu et popSimu.select)
affichDataH <- function(tirages,Nb_rep){

  df <- NULL
  for(i in 1:Nb_rep)
  {
    len <- tirages[[i]]$fixationTime
    temp_df <- data.frame(x=0:len, y=tirages[[i]]$coeffH, col=rep(i, each=len+1))
    df <- rbind(df,temp_df)
  }
  ggplot(df,aes(x=x,y=y,group=col,colour=factor(col))) + geom_line(size=1) + labs(title="Value of Heterozysity ratio over Time",x="Time", y = "CoeffH") +
    theme(axis.title=element_text(size=14,face="bold"),plot.title=element_text(size=18,face="bold",hjust=0.5),legend.position = "none") # plot data

}

#histogramme du temps de fixation (pour popSimu et popSimu.select)
affichDataT <- function(vect_temps,bw){
  colors <- c("mean" = "blue", "variance" = "red")
  texte <- c(paste("Mean = ",round(mean(vect_temps),2)), paste("Variance = ", round(var(vect_temps),2)))
  df <- as.data.frame(vect_temps)
  ggplot(df,aes(x=vect_temps)) + geom_histogram(color="blue",fill="lightblue",binwidth=bw)+
    geom_vline(aes(xintercept=mean(vect_temps),color="mean"), linetype="dashed", size=1)+
    geom_vline(aes(xintercept=var(vect_temps),color="variance"), linetype="dashed", size=1)+
    labs(title="Fixation time histogram plot",x="Fixation time", y = "Count")+
    theme(axis.title=element_text(size=14,face="bold"),plot.title=element_text(size=18,face="bold",hjust=0.5) )+
    scale_color_manual(name = "Légende", labels=texte, values = colors)+
    coord_cartesian(xlim = c(min(vect_temps), max(vect_temps)))
}

################### temps de fixation en fonction des probas initiales ###################

## pour simu.simu  (plusieurs simulations de popSimu)
#temps de fixation sans sélection + cas discret + différences
affichTP <- function (listSimu,ne,step){

  vecteurTemps <- sapply((lapply(listSimu, vectFix)), mean) #pour chaque mclapply, moyenner les Temps de fixation obtenus
  proba <- vectPropa(listSimu)

  df <- data.frame(proba=vectPropa(listSimu),time=vecteurTemps)

  ## courbe théorique ##
  FT <- -4*ne*(proba*log(proba) + (1-proba)*log(1-proba))
  FT[1] <- 0
  FT[length(FT)] <- 0
  vect_proba <- seq(0,1, by=0.001)
  FT_vect_proba <- -4*ne*(vect_proba*log(vect_proba) + (1-vect_proba)*log(1-vect_proba))
  FT_vect_proba[1] <- 0
  FT_vect_proba[length(FT_vect_proba)] <- 0
  df2 <- data.frame(proba= vect_proba, time=FT_vect_proba) #df des proba au 1/100

  ###### cas discret #####
  mat_un <- rep(1,(2*ne)-1)
  Im <- diag(1, (2*ne)-1) #matrice identité

  p_ij <- NULL #vecteur des probas de passages de i à j

  for (i in 1:(2*ne-1)){
    for (j in 1:(2*ne-1)){
      p_ij <- c(p_ij,dbinom(x=j,size=2*ne,prob=i/(2*ne)))
    }
  }

  mat_pij <- matrix(p_ij,(2*ne)-1,byrow=TRUE) #remplir la matrice en ligne
  m <- rowSums(solve(Im - mat_pij)) #déterminer m
  m <- c(0,m,0) #compléter avec les valeurs manquantes

  proba.m <- (0:(2*ne))/(2*ne)
  df3 <- data.frame(proba = proba.m, time=m) #discret

  temps.theo.discret <- -4*ne*(proba.m*log(proba.m) + (1-proba.m)*log(1-proba.m))
  temps.theo.discret[1] <- 0
  temps.theo.discret[length(temps.theo.discret)] <- 0
  df4 <- data.frame(proba=proba.m, time=temps.theo.discret)

  ##affichage
  colors <- c("simulation" = "blue", "théorie" = "red", "discret" = "#C2F732", "simu-theo" = "#FF00FF", "discret-theo"="yellow")
  ggplot(df,aes(x=proba,y=time,color="simulation")) + geom_point() +
    labs(title="Fixation time in terms of inital probability",x="p", y = "Fixation time")+
    theme(axis.title=element_text(size=14,face="bold"),plot.title=element_text(size=18,face="bold",hjust=0.5)) +
    geom_smooth(linetype="blank", level=0.95,size=0.8,method = 'loess', formula= y ~ x) +
    geom_line(data=df2,aes(x=proba, y=FT_vect_proba,color="théorie"), size=1.2,) +
    geom_line(aes(x=proba, y=(time-FT),color="simu-theo"), size=1.2) +
    geom_point(data=df3, aes(x=proba,y=time, color="discret"), size=1.2) +
    geom_line(data=df4, aes(x=proba, y=temps.theo.discret-m,color="discret-theo"), size=1.2)+
    scale_color_manual(name = "Légende", values = colors)
}

## pour simu.simu.select (plusieurs simulations de popSimu.select)
## non utilisé car on utilise affichTP.both
affichTP.select <- function (listSimu){
  vecteurTemps <- sapply((lapply(listSimu, vectFix)), mean) #pour chaque mclapply, moyenner les Temps de fixation obtenus
  df <- data.frame(proba=vectPropa(listSimu),time=vecteurTemps)

  ggplot(df,aes(x=proba,y=time)) + geom_point(color="blue", size=1)  +
    labs(title="Fixation time in terms of inital probability",x="p", y = "Fixation time")+
    theme(axis.title=element_text(size=14,face="bold"),plot.title=element_text(size=18,face="bold",hjust=0.5)) +
    geom_smooth(linetype="blank", level=0.95,size=0.8,method = 'loess', formula= y ~ x) +
    geom_line(aes(x=proba, y=time),color="blue")
}

## pour simu.simu.select (plusieurs simulations de popSimu.select)
##temps de fixation avec et sans sélection + cas discret + différences
affichTP.both <- function (listSimu,listSimuS,ne,s){

  f <- function(x)
  {return(exp(-2*alpha*x)/(2*alpha))}

  g <- function(x)
  {return((2*ne/(alpha))*(exp(-2*alpha*x)*expint_Ei(2*alpha*x) - exp(-2*alpha*(x-1))* expint_Ei(2*alpha*(x-1)) + log(1/x-1)))}

  #théorie avec sélection
  colors <- c("sélection"="blue", "sans sélection"="black", "discret"="#16B84E", "théorie"="red", "discret-théorie" = "yellow", "simu-théorie" = "#FF00FF", "théorie 2" = "purple")
  alpha <- s*2*ne
  step <- 0.001
  pTheo <- seq(step,1-step,by =step)
  gamma <- 0.577215664901532
  k <-2
  C <- (g(1-10^-10)- g(10^-10))/(f(1)-f(0))
  K <- C*f(1) -g(1- 10^-10)

  m_p <- K - C*exp(-2*alpha*pTheo)/(2*alpha) + (2*ne/alpha)*(exp(-2*alpha*pTheo)*expint_Ei(2*alpha*pTheo) - exp(-2*alpha*(pTheo-1))* expint_Ei(2*alpha*(pTheo-1)) + log(1/pTheo-1))
  #m_pApprox <-  K - C*exp(-2*alpha*pTheo)/(2*alpha) + (2*ne/alpha) * (exp(-2*alpha*pTheo)* (1-exp(2*alpha)) * (gamma - integrale2(0.01, 2*alpha*(pTheo-1),k)) - exp(-2*alpha*pTheo) * integrale2(2*alpha*pTheo,2*alpha*(pTheo-1), k))
  #m_pApprox <-  K - C*exp(-2*alpha*pTheo)/(2*alpha) + (2*ne/alpha) * (-exp(-2*alpha*pTheo)* integrale2(0.01, 2*alpha*pTheo,k)- exp(2*alpha*(1-pTheo)) * integraleNeg(0.01, 2*alpha*(1-pTheo),k) + gamma * exp(-2*alpha*pTheo)* (1-exp(2*alpha)))
  #m_pApprox <- K - C*exp(-2*alpha*pTheo)/(2*alpha) + (2*ne/alpha)*(- exp(-2*alpha*p) * integrale2(0.01, 2*alpha*p, k) + exp(-2*alpha*p)*gamma + log(2*alpha))

  pTheo <- c(0, pTheo, 1)
  m_p <- c(0, m_p, 0)

  #m_pApprox <- c(0,m_pApprox,0)
  df.theo <- data.frame(proba=pTheo, temps=m_p)
  #df.theoApprox <- data.frame(proba=pTheo, temps=m_pApprox)

  #sans sélection
  #vecteurTemps <- sapply((lapply(listSimu, FUN=vectFix)), FUN=mean)
  #df <- data.frame(proba=vectPropa(listSimu),time=vecteurTemps)

  #théorie sans sélection
  vect_proba <- seq(0,1, by=0.001)
  FT_vect_proba <- -4*ne*(vect_proba*log(vect_proba) + (1-vect_proba)*log(1-vect_proba))
  FT_vect_proba[1] <- 0
  FT_vect_proba[length(FT_vect_proba)] <- 0
  df.theo.noSelect <- data.frame(proba= vect_proba, time=FT_vect_proba) #df des proba au 1/1000


  #avec sélection
  vecteurTempsS <- sapply((lapply(listSimuS, FUN=vectFix)), FUN=mean)
  dfS <- data.frame(probaS=vectPropa(listSimuS),timeS=vecteurTempsS)

  ###### cas discret #####
  mat_un <- rep(1,(2*ne)-1)
  Im <- diag(1, (2*ne)-1) #matrice identité

  p_ij <- NULL #vecteur des probas de passages de i à j

  for (i in 1:(2*ne-1)){
    for (j in 1:(2*ne-1)){
      p_ij <- c(p_ij,dbinom(x=j,size=2*ne,prob= (1+s)*i/((1+s)*i + 2*ne - i)))
    }
  }

  mat_pij <- matrix(p_ij,(2*ne)-1,byrow=TRUE) #remplir la matrice en ligne
  m <- rowSums(solve(Im - mat_pij)) #déterminer m
  m <- c(0,m,0) #compléter avec les valeurs manquantes
  proba.m <- (0:(2*ne))/(2*ne)
  df3 <- data.frame(proba = proba.m, time=m) #discret

  temps.theo.discret <- -4*ne*(proba.m*log(proba.m) + (1-proba.m)*log(1-proba.m))
  temps.theo.discret[1] <- 0
  temps.theo.discret[length(temps.theo.discret)] <- 0
  df4 <- data.frame(proba=proba.m, time=temps.theo.discret)
  ###### fin cas discret #####

  ## discret - théorie ##
  proba.m <- proba.m[2:(length(proba.m)-1)]
  m_p2 <-  K - C*exp(-2*alpha*proba.m)/(2*alpha) + (2*ne/alpha)*(exp(-2*alpha*proba.m)*expint_Ei(2*alpha*proba.m) - exp(-2*alpha*(proba.m-1))* expint_Ei(2*alpha*(proba.m-1)) + log(1/proba.m-1))
  proba.m <- c(0,proba.m,1)
  m_p2 <- c(0,m_p2,0)
  df.discret.theo <- data.frame(proba=proba.m, time= m-m_p2 )

  ## simulation - théorie ##
  vectProbaListSimuS <- vectPropa(listSimuS)[2:(length(vectPropa(listSimuS))-1)]
  temps.simu.theo <- K - C*exp(-2*alpha*vectProbaListSimuS)/(2*alpha) + (2*ne/alpha)*(exp(-2*alpha*vectProbaListSimuS)*expint_Ei(2*alpha*vectProbaListSimuS) - exp(-2*alpha*(vectProbaListSimuS-1))* expint_Ei(2*alpha*(vectProbaListSimuS-1)) + log(1/vectProbaListSimuS-1))
  temps.simu.theo <- c(0,temps.simu.theo,0)
  df.simu.theo <- data.frame(proba=vectPropa(listSimuS), time= vecteurTempsS - temps.simu.theo)

  ggplot(dfS,aes(x=probaS, y=timeS,color="sélection")) +  geom_point() +
    labs(title="Fixation time in terms of inital probability",x="p", y = "Fixation time")+
    theme(axis.title=element_text(size=14,face="bold"),plot.title=element_text(size=18,face="bold",hjust=0.5)) +
    geom_line(data=df.theo.noSelect, aes(x=proba, y=time, colour="sans sélection"),size=1) +
    geom_point(data=df3, aes(x=proba,y=time, color="discret"), size=1.4) +
    geom_line(data=df.theo, aes(x=proba, y=temps, color="théorie")) +
    #geom_line(data=df.theoApprox, aes(x=proba, y=temps, color="théorie 2")) +
    geom_line(data= df.discret.theo, aes(x=proba, y=time, color="discret-théorie")) +
    geom_line(data= df.simu.theo, aes(x=proba, y=time, color="simu-théorie")) +
    scale_color_manual(name = "Légende", values = colors)

  #ggplot(df.theoApprox, aes(x=proba, y=temps)) + geom_line()
}

##variance du temps de fixation
affichVP <- function (listSimu){
  colors <- c("variance"="blue", "courbe de tendance"="#FFFF6B")
  vecteurTemps <- sapply((lapply(listSimu, vectFix)), var) #pour chaque mclapply, retourner la variance les Temps de fixation obtenus
  df <- data.frame(proba=vectPropa(listSimu),time=vecteurTemps)
  ggplot(df,aes(x=proba,y=time,color="variance")) + geom_line()  + labs(title=" Variance in terms of inital probability",x="p", y = "Variance")+
    theme(axis.title=element_text(size=14,face="bold"),plot.title=element_text(size=18,face="bold",hjust=0.5)) +
    geom_smooth(aes(color="courbe de tendance"), level=0.95,size=0.8,method = 'loess', formula= y ~ x)+
    scale_color_manual(name = "Légende", values = colors)

}

###########

## c'est comme la fonction mclapply() de R mais faite à la main
## on ne s'en sert pas

mclapply.simu <- function(Ne,Na, nbRep){
  listout <- vector("list", nbRep)
  for(i in 1:nbRep){
    listout[[i]] <- popSimu(Ne,Na)
  }
  return(listout)
}

mclapply.simu.select <- function(Ne,Na, nbRep,s){
  listout <- vector("list", nbRep)
  for(i in 1:nbRep){
    listout[[i]] <- popSimu.select(Ne,Na,s)
  }
  return(listout)
}

##############################################
## plusieurs popSimu avec et sans sélection pour toutes les proba initiales calculées en fonction de 2*Ne et step

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
  if(Ne%%step != 0){
    valA <- c(valA,2*Ne)
  }
  print(valA)
  outF <- vector("list", length(valA))
  for (i in 1:length(valA)){
    outF[[i]] <- mclapply(rep(Ne,Nb_rep),popSimu.select, nbA = valA[i],s=s,mc.cores = 1)
  }
  return(outF)
}
