
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

