plotH <- function(generation){
  plot(generation$coeffH, type="l", ylim= c(min(generation$coeffH), max(generation$coeffH)),
       xlab="Temps de fixation", ylab="Coeff d'hétérozygotie")
}

plotA <- function(generation){
  plot(generation$nbA, type="l", xlim=c(0,generation$fixationTime),ylim= c(min(generation$nbA),
                                                                           max(generation$nbA)),
       xlab="Temps de fixation", ylab="Nombre d'allèles A")
  abline(h=c(min(generation$nbA),max(generation$nbA)), col=2)
}

affichDataA <- function(tirages,Nb_rep){

  df <- NULL
  for(i in 1:Nb_rep)
  {
    len <- tirages[[i]]$fixationTime
    temp_df <- data.frame(x=1:len, y=tirages[[i]]$nbA, col=rep(i, each=len))
    df <- rbind(df,temp_df)
  }
  ggplot(df,aes(x=x,y=y,group=col,colour=factor(col))) + geom_line() # plot data

}


affichDataH <- function(tirages,Nb_rep){

  df <- NULL
  for(i in 1:Nb_rep)
  {
    len <- tirages[[i]]$fixationTime
    #print(len)
    temp_df <- data.frame(x=1:len, y=tirages[[i]]$coeffH, col=rep(i, each=len))
    df <- rbind(df,temp_df)
  }
  ggplot(df,aes(x=x,y=y,group=col,colour=factor(col))) + geom_line() # plot data

}

affichDataT <- function(vect_temps,bw){
  df <- as.data.frame(vect_temps)
  ggplot(df,aes(x=vect_temps)) + geom_histogram(color="blue",fill="lightblue",binwidth=bw)+ geom_vline(aes(xintercept=mean(vect_temps)),color="blue", linetype="dashed", size=1)+
    labs(title="Fixation time histogram plot",x="Fixation time", y = "Count")+theme(axis.title=element_text(size=14,face="bold"),plot.title=element_text(size=18,face="bold",hjust=0.5))
}
