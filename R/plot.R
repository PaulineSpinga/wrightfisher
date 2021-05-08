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
