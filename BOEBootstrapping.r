#rm(list=ls())
#install.packages('MCMCpack')
library(MCMCpack)
library(NSM3)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Rubin's Bayesian bootstrap for mean difference
mean.difference.rbb<-function(y1,y2,n)
{
  
  y1<-sort(y1)
  y2<-sort(y2)  
  
  draw<-matrix(nrow=n,ncol=1,0)
 
  for (i in 1:n)
  {
    
    u1<-c(0,runif(n=length(y1)-1),1)
    u1<-u1[order(u1)]
    g1<-u1
    
    u2<-c(0,runif(n=length(y2)-1),1)
    u2<-u2[order(u2)]    
    g2<-u2
    
    for(j in 2:length(u1)-1)
    {
      g1[j]<-u1[j+1]-u1[j]
    }

    for(k in 2:length(u2)-1)
    {
      g2[k]<-u2[k+1]-u2[k]
    }    

    g1<-g1[2:length(g1)-1]    
    g2<-g2[2:length(g2)-1]
    
    draw[i,1]<-weighted.mean(x=y1,w=g1,na.rm=TRUE)-weighted.mean(x=y2,w=g2,na.rm=TRUE)
    
  }
  return(draw)
}


rubin.dirichlet<-function(n,bins)
{
  g1<-matrix(nrow=n,ncol=bins)
  
  for(i in 1:n)
  {
    u1<-c(0,runif(n=bins-1),1)
    u1<-u1[order(u1)]
  

      for(j in 2:length(u1)-1)
      {
        g1[i,j]<-u1[j+1]-u1[j]
      }
    }
return(g1) 
}


bb.simple <- function(y,n)
{
  bb<-apply(rdirichlet(n,rep(1,length(y))),1,weighted.mean,x=y,na.rm=T)
  return(bb)
}

bb.DiD<- function(y1,y2,n)
{
  bb<-apply(rdirichlet(n,rep(1,length(y1))),1,weighted.mean,x=y1,na.rm=T)-apply(rdirichlet(n,rep(1,length(y2))),1,weighted.mean,x=y2,na.rm=T)
  return(bb)
}

bb.DiD<- function(y1,y2,n)
{
  bb<-apply(rdirichlet(n,rep(1,length(y1))),1,weighted.mean,x=y1,na.rm=T)-apply(rdirichlet(n,rep(1,length(y2))),1,weighted.mean,x=y2,na.rm=T)
  return(bb)
}

#ferg.df
#function (x, alpha, mu, npoints, ...) 
#{
#  x = sort(x)
#  n = length(x)
#  Fn = ecdf(x)
#  temp.x = seq(min(x), max(x), length.out = npoints)
#  f.n = Fn(temp.x)
#  f.0 = mu(temp.x, ...)
#  mu.n = (alpha/(alpha + n)) * f.0 + (n/(alpha + n)) * f.n
#  return(mu.n)
#}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#install.packages('ecdf')

#file.choose()

forestPlotData<-data.frame(read.csv( "C:\\Users\\mholtman\\Documents\\GitHub\\BackOfTheEnvelopeImpact\\forestPlotData.csv",fileEncoding="UTF-8-BOM",stringsAsFactors=F))
data<-forestPlotData$LogOdds

marginal.change<-function(x)
  {
  output<-matrix(nrow=length(x),ncol=1,0)
  for(i in 2:length(x))
  {
    output[i-1]<-(x[i]-x[i-1])
  }
  return(output)
 }

mu<-function(x){return(pnorm(sort(x),mean=0,sd=2,lower.tail=T,log.p=F))}

replicates<-bb.simple(forestPlotData$LogOdds,n=1000)
update<-ferg.df(x=sort(replicates),alpha=1,mu=mu,npoints=10)

x.even<-seq(min(replicates), max(replicates), length.out = 10)
tick.gap<-(x.even[length(x.even)]-x.even[length(x.even)-1])/2
x.even.expanded<-c(x.even,x.even[length(x.even)]+tick.gap)

ferg.df.marginal<-data.frame(x.even,marginal.change(update))

h<-hist(replicates,breaks=x.even.expanded)
counts<-data.frame(h$mids,h$density)
colnames(counts)<-c('value','prior')

data<-rnorm(mean=0,sd=2,n=1000)
#h<-hist(data,breaks=x.even.expanded,freq=F)
h<-hist(data[data>=min(x.even.expanded)&data<=max(x.even.expanded)],breaks=x.even.expanded,freq=F)
#h<-hist(data,breaks=seq(x.even.expanded[length(x.even.expanded)],min(data),-tick.gap),freq=F)
counts2<-data.frame(h$mids,h$density)
colnames(counts2)<-c('value','data')
combined<-merge(counts,counts2)
combined

n<-1000
alpha<-1000
combined$posterior<-(n/(alpha+n))*(combined$prior/sum(combined$prior))+(alpha/(alpha+n))*(combined$data/sum(combined$data))

combined$tick<-combined$value-tick.gap
plot(combined$tick,combined$posterior,type='l')
points(ferg.df.marginal$x.even,ferg.df.marginal$marginal.change.update,type='l',col='red')

#hist(replicates,xlim=c(-4,4),breaks=20)
#hist(data)

plot(density(replicates))
points(density(data),type='l',col='blue')
points(combined$tick,combined$posterior,type='l',col='red')
points(ferg.df.marginal$x.even,ferg.df.marginal$marginal.change.update,type='l',col='red')

plot(density(data))

combined
ferg.df.marginal
