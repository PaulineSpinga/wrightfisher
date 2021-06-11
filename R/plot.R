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

affichDataT <- function(vect_temps,bw){
  colors <- c("mean" = "blue", "variance" = "red")
  texte <- c(paste("Mean = ",round(mean(vect_temps),2)), paste("Variance = ", round(var(vect_temps),2)))
  print(paste("Variance = ", round(var(vect_temps),2)))
  df <- as.data.frame(vect_temps)
  ggplot(df,aes(x=vect_temps)) + geom_histogram(color="blue",fill="lightblue",binwidth=bw)+
    geom_vline(aes(xintercept=mean(vect_temps),color="mean"), linetype="dashed", size=1)+
    geom_vline(aes(xintercept=var(vect_temps),color="variance"), linetype="dashed", size=1)+
    labs(title="Fixation time histogram plot",x="Fixation time", y = "Count")+
    theme(axis.title=element_text(size=14,face="bold"),plot.title=element_text(size=18,face="bold",hjust=0.5) )+
    scale_color_manual(name = "Légende", labels=texte, values = colors)

}


affichTP <- function (listSimu,ne,step){

  vecteurTemps <- sapply((lapply(listSimu, vectFix)), mean) #pour chaque mclapply, moyenner les Temps de fixation obtenus
  proba <- vectPropa(listSimu)

  df <- data.frame(proba=vectPropa(listSimu),time=vecteurTemps)

  FT <- -4*ne*(proba*log(proba) + (1-proba)*log(1-proba))
  FT[1] <- 0
  FT[length(FT)] <- 0
  vect_proba <- seq(0,1, by=0.001)
  FT_vect_proba <- -4*ne*(vect_proba*log(vect_proba) + (1-vect_proba)*log(1-vect_proba))
  FT_vect_proba[1] <- 0
  FT_vect_proba[length(FT_vect_proba)] <- 0
  df2 <- data.frame(proba= vect_proba, time=FT_vect_proba) #df des proba au 1/100

  mat_un <- rep(1,(2*ne)-1)
  Im <- diag(1, (2*ne)-1)

  p_ij <- NULL

  for (i in 1:(2*ne-1)){
    for (j in 1:(2*ne-1)){
      p_ij <- c(p_ij,dbinom(x=j,size=2*ne,prob=i/(2*ne)))
    }
  }

  mat_pij <- matrix(p_ij,(2*ne)-1,byrow=TRUE)
  m <- rowSums(solve(Im - mat_pij))
  m <- c(0,m,0)

  proba.m <- (0:(2*ne))/(2*ne)
  df3 <- data.frame(proba = proba.m, time=m)

  temps.theo.discret <- -4*ne*(proba.m*log(proba.m) + (1-proba.m)*log(1-proba.m))
  temps.theo.discret[1] <- 0
  temps.theo.discret[length(temps.theo.discret)] <- 0
  df4 <- data.frame(proba=proba.m, time=temps.theo.discret)

  ##affichage
  colors <- c("simulation" = "blue", "théorie" = "red", "discret" = "#C2F732", "simu-theo" = "#FF00FF", "discret-theo"="yellow")
  ggplot(df,aes(x=proba,y=time,color="simulation")) + geom_point() +
    labs(title="Fixation time in terms of inital probability",x="p", y = "Fixation time")+
    theme(axis.title=element_text(size=14,face="bold"),plot.title=element_text(size=18,face="bold",hjust=0.5)) +
    geom_smooth(linetype="blank", level=0.95,size=0.8) +
    geom_line(data=df2,aes(x=proba, y=FT_vect_proba,color="théorie"), size=1.2,) +
    geom_line(aes(x=proba, y=(time-FT),color="simu-theo"), size=1.2) +
    geom_point(data=df3, aes(x=proba,y=time, color="discret"), size=1.2) +
    geom_line(data=df4, aes(x=proba, y=temps.theo.discret-m,color="discret-theo"), size=1.2)+
    scale_color_manual(name = "Légende", values = colors)
}



affichTP.select <- function (listSimu){
  vecteurTemps <- sapply((lapply(listSimu, vectFix)), mean) #pour chaque mclapply, moyenner les Temps de fixation obtenus
  df <- data.frame(proba=vectPropa(listSimu),time=vecteurTemps)

  ggplot(df,aes(x=proba,y=time)) + geom_point(color="blue", size=1)  +
    labs(title="Fixation time in terms of inital probability",x="p", y = "Fixation time")+
    theme(axis.title=element_text(size=14,face="bold"),plot.title=element_text(size=18,face="bold",hjust=0.5)) +
    geom_smooth(linetype="blank", level=0.95,size=0.8) +
    geom_line(aes(x=proba, y=time),color="blue")
}

affichTP.both <- function (listSimu,listSimuS,ne,s){

  f <- function(x)
  {return(exp(-2*alpha*x)/(2*alpha))}

  g <- function(x)
  {return((2*ne/(alpha))*(exp(-2*alpha*x)*expint_Ei(2*alpha*x) - exp(-2*alpha*(x-1))* expint_Ei(2*alpha*(x-1)) + log(1/x-1)))}

  colors <- c("sélection"="blue", "sans sélection"="black", "discret"="#16B84E", "théorie"="red", "discret-théorie" = "yellow", "simu-théorie" = "#FF00FF")
  alpha <- s*2*ne
  step <- 0.001
  pTheo <- seq(step,1-step,by =step)
  gamma <- 0.577215664901532
  C <- (g(1-10^-10)- g(10^-10))/(f(1)-f(0))
  K <- C*f(1) -g(1- 10^-10)

  m_p <- K - C*exp(-2*alpha*pTheo)/(2*alpha) + (2*ne/alpha)*(exp(-2*alpha*pTheo)*expint_Ei(2*alpha*pTheo) - exp(-2*alpha*(pTheo-1))* expint_Ei(2*alpha*(pTheo-1)) + log(1/pTheo-1))

  pTheo <- c(0, pTheo, 1)
  m_p <- c(0, m_p, 0)
  df.theo <- data.frame(proba=pTheo, temps=m_p)


  vect_proba <- seq(0,1, by=0.001)
  FT_vect_proba <- -4*ne*(vect_proba*log(vect_proba) + (1-vect_proba)*log(1-vect_proba))
  FT_vect_proba[1] <- 0
  FT_vect_proba[length(FT_vect_proba)] <- 0
  df.theo.noSelect <- data.frame(proba= vect_proba, time=FT_vect_proba) #df des proba au 1/1000

  vecteurTempsS <- sapply((lapply(listSimuS, FUN=vectFix)), FUN=mean)
  dfS <- data.frame(probaS=vectPropa(listSimuS),timeS=vecteurTempsS)

  mat_un <- rep(1,(2*ne)-1)
  Im <- diag(1, (2*ne)-1)

  p_ij <- NULL

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


  m_p2 <-  m_p <- K - C*exp(-2*alpha*proba.m)/(2*alpha) + (2*ne/alpha)*(exp(-2*alpha*proba.m)*expint_Ei(2*alpha*proba.m) - exp(-2*alpha*(proba.m-1))* expint_Ei(2*alpha*(proba.m-1)) + log(1/proba.m-1))
  m_p2[1] <- 0
  m_p2[length(m_p2)] <- 0
  df.discret.theo <- data.frame(proba=proba.m, time= m-m_p2 )

  ## simulation - théorie ##
  temps.simu.theo <- m_p <- K - C*exp(-2*alpha*vectPropa(listSimuS))/(2*alpha) + (2*ne/alpha)*(exp(-2*alpha*vectPropa(listSimuS))*expint_Ei(2*alpha*vectPropa(listSimuS)) - exp(-2*alpha*(vectPropa(listSimuS)-1))* expint_Ei(2*alpha*(vectPropa(listSimuS)-1)) + log(1/vectPropa(listSimuS)-1))
  temps.simu.theo[1] <- 0
  temps.simu.theo[length(temps.simu.theo)] <- 0
  df.simu.theo <- data.frame(proba=vectPropa(listSimuS), time= vecteurTempsS - temps.simu.theo)


  ggplot(dfS,aes(x=probaS, y=timeS,color="sélection")) +  geom_point() +
    labs(title="Fixation time in terms of inital probability",x="p", y = "Fixation time")+
    theme(axis.title=element_text(size=14,face="bold"),plot.title=element_text(size=18,face="bold",hjust=0.5)) +
    geom_line(data=df.theo.noSelect, aes(x=proba, y=time, colour="sans sélection"),linetype="dotted",size=1.2) +
    geom_point(data=df3, aes(x=proba,y=time, color="discret"), size=1.4) +
    geom_line(data=df.theo, aes(x=proba, y=temps, color="théorie")) +
    geom_line(data= df.discret.theo, aes(x=proba, y=time, color="discret-théorie")) +
    geom_line(data= df.simu.theo, aes(x=proba, y=time, color="simu-théorie")) +
    scale_color_manual(name = "Légende", values = colors)
}


affichVP <- function (listSimu){
  colors <- c("variance"="blue", "courbe de tendance"="#FFFF6B")
  vecteurTemps <- sapply((lapply(listSimu, vectFix)), var) #pour chaque mclapply, retourner la variance les Temps de fixation obtenus
  df <- data.frame(proba=vectPropa(listSimu),time=vecteurTemps)
  ggplot(df,aes(x=proba,y=time,color="variance")) + geom_line()  + labs(title=" Variance in terms of inital probability",x="p", y = "Variance")+
    theme(axis.title=element_text(size=14,face="bold"),plot.title=element_text(size=18,face="bold",hjust=0.5)) +
    geom_smooth(aes(color="courbe de tendance"), level=0.95,size=0.8)+
    scale_color_manual(name = "Légende", values = colors)

}
