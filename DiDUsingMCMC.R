library(rjags)

#file.choose()
experiment<-data.frame(read.csv("C:\\Users\\mchol\\Documents\\Writing\\Bayes\\Data\\buitrago.csv",stringsAsFactors = F))
experiment$difference<-experiment$RelativeRisk
experiment$arm<-0
experiment$arm[experiment$Outcome == 'Any cardiovascular disease']<-1

set.seed(100)
n<-100
x<-rnorm(n,0,5)


model.string<-
"
model
{

  #Data
  for (i in 1:N)
  {
    Y[i]~dnorm(yhat[i],tau)
    yhat[i]=a+b*treatment[i]
  }
  
  #Priors
  a~dunif(-2,2)
  b~dunif(-2,2)

  tau<-pow(sigma,-2)
  sigma~dunif(0,100)
}
" 

model.spec<-textConnection(model.string)

full.model<-jags.model(file=model.spec
                 ,data=list( 'Y'=experiment$difference
                            ,'treatment'=experiment$arm
                            ,'N'=nrow(experiment)
                            ),n.chains=4,n.adapt=100)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

update(full.model,1000)

jags.samples(full.model,c('b','tau'),1000)

#Now sample from the posterior
samples<-coda.samples(full.model,variable.names = c('a','b'),n.iter = 10000)
summary(samples)
plot(samples)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

simpleModel<-lm(data=experiment,difference~arm)
summary(simpleModel)
