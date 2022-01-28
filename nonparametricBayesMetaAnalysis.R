library(NSM3)
library(MCMCpack)
#file.choose()
#install.packages('MCMCpack')

#The multinomial distribution takes a vector of probabilities and returns counts:
p<-c(0.15,.2,.55,.1)
#Make sure the probabilities form a proper distribution: 
sum(p)
rmultinom(n=1,p=p,size=100)

#Whereas the dirichlet distribution takes counts and produces probabilities.
inputs<-rmultinom(n=1,p=c(.33,.33,.33),size=10)
rdirichlet(n=10,alpha=inputs)

#The improper prior D(0,0,0) doesn't give us anything.
#Gelman says it's uniform in the log (thetas). What does that mean?
rdirichlet(n=1,alpha=c(0,0,0))

sum(inputs)


#trials<-data.frame(read.csv("C:\\Users\\mchol\\Documents\\Writing\\Bayes\\TrialData.csv"))
trials<-data.frame(read.csv("C:\\Users\\mchol\\Documents\\Writing\\Bayes\\data\\hofmann.csv"))

trials
#x<-trials$Difference

backbone<-function(x,npoints){return(seq(min(x), max(x), length.out = npoints))}

cunif<-function(x)
{
  data<-sort(x)
  return((data-min(data))/(max(data)-min(data)))
}

plot(sort(trials$EffectSize),cunif(trials$EffectSize),type='l')

prior<-cunif(trials$EffectSize)

diff(prior)

plot(sort(trials$EffectSize)[2:nrow(trials)],diff(prior),type='l')



bayesBB<-function(x,n)
{
  output<-matrix(nrow=n,ncol=1)
  
  for(i in 1:n)
  {
    output[i,1]<-weighted.mean(x,w=rdirichlet(n=x,alpha=rep(1,length(x))))
  }
  return(output)
}

BBFunction<-function(x,f,n)
{
  output<-matrix(nrow=n,ncol=1)
  
  for(i in 1:n)
  {
    weighted.x<-x*rdirichlet(n=x,alpha=rep(1,length(x)))
    output[i,1]<-f(weighted.x)
  }
  return(output)
}


data<-trials$EffectSize

hist(BBFunction(x=data,f=mean,n=1000),breaks=50)
hist(BBFunction(x=data,f=sd,n=1000),breaks=50)


#getAnywhere(ecdf)

# ecdf<-function (x) 
# {
#   x <- sort(x)
#   n <- length(x)
#   if (n < 1) 
#     stop("'x' must have 1 or more non-missing values")
#   vals <- unique(x)
#   rval <- approxfun(vals, cumsum(tabulate(match(x, vals)))/n, 
#                     method = "constant", yleft = 0, yright = 1, f = 0, ties = "ordered")
#   class(rval) <- c("ecdf", "stepfun", class(rval))
#   assign("nobs", n, envir = environment(rval))
#   attr(rval, "call") <- sys.call()
#   rval
# }
# 
# 
# ferg.df<-function (x, alpha, mu, npoints, ...) 
# {
#   x = sort(x)
#   n = length(x)
#   Fn = ecdf(x)
#   temp.x = seq(min(x), max(x), length.out = npoints)
#   f.n = Fn(temp.x)
#   f.0 = mu(temp.x, ...)
#   mu.n = (alpha/(alpha + n)) * f.0 + (n/(alpha + n)) * f.n
#   return(mu.n)
# }

#observed<-trials$Difference
#observed<-BBBDiff

dev.off()

npoints<-100

trialData<-sort(trials$EffectSize)

#cutoffs<-backbone(trials$Difference,npoints)
cutoffs<-backbone(trials$EffectSize,npoints)

Fn = ecdf(trials$EffectSize)
observed<-Fn(cutoffs)

#prior<-cunif(trialData,min(trialData),max(trialData))
prior<-cunif(cutoffs)

BBDraws<-bayesBB(trials$EffectSize,1000)
Fn = ecdf(BBDraws)
BBPosterior<-Fn(cutoffs)
hist(BBDraws,breaks=50)

length(trialData)
posterior<-ferg.df(x=trialData,alpha=35,npoints=npoints,mu=dunif,min(cutoffs),max(cutoffs))

output<-data.frame(cutoffs,prior,observed,posterior,BBPosterior)
plot(output$cutoffs,output$prior,type='l',col='yellow',xlab='Effect size',ylab='Probability')
lines(output$cutoffs,output$observed,type='l',col='red')
lines(output$cutoffs,output$posterior,type='l',col='black')

output
max(output$cutoffs[which(output$prior<=.5)])
max(output$cutoffs[which(output$posterior<=.5)])


