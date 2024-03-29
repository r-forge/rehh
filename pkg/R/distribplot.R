distribplot<-function(data,col=c("blue","red"),main="iHS distribution",xlab="iHS"){
     plot.density(density(data,na.rm=T),main=main,xlab="iHS",col=col[1])
     curve(dnorm,col=col[2],add=T)
     smartlegend("right","top",c("Observed","Gaussian"),fill=col)
  }
