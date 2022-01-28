install.packages('BayesFactor')
require(rjags)
library(rjags)
library(BayesFactor)

set.seed(100)
n<-100
x<-rnorm(n,0,5)

model.string<-
"
model 
{

  for (i in 1:N){
  x[i]~dnorm(mu,tau)
  }
  
  mu ~ dnorm(0,.0001)
  tau<- pow(sigma,-2)
  sigma ~ dunif(0,100)
}
" 

model.spec<-textConnection(model.string)



jags<-jags.model(file=model.spec,data=list('x'=x,'N'=n),n.chains=4,n.adapt=100)


update(jags,1000)


jags.samples(jags,c('mu','tau'),1000)
