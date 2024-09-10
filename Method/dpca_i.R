#Get FPC scores for single sample in single sensor 
#X2 has 50 time points, but sample# = 1, and p = 1
dpca_score<-function(X2,eval,efns,mua){ 
  demeaned <- X2-mua
  n <- dim(efns)[1] 
  d = length(eval)
  scores <-rep(0,d)
  for(l in 1:d){
    scores[l] = sum(demeaned*efns[,l])/n
  }
  scores #return a d*1 matrix
}
#Estimate fpca parameters 
#input x have samples and time points but only 1 sensor
#d0 is the number of decomposition dimension.
dpca_est<-function(x,d0){
  m = dim(x)[1]
  n = dim(x)[2]
  cova<-array(0,dim=c(n,n))
  mean.hat<-array(0,dim=c(n,1))
  mean.hat <- colMeans(x)
  #construct covariance fucntion 
  for (i in 1:m){
    cova = cova + (x[i,]-mean.hat)%*%t(x[i,]-mean.hat)
  }
  cova = (1/m)*cova
  #Eigen decompostion 
  eio<-eigen(cova)
  eval<-eio$values/n
  evec<-eio$vectors*sqrt(n)
  #truncate for d = d0
  eval0<-eval[1:d0]  
  efns0<-evec[,1:d0]
  demeaned <- x - t(matrix(rep(mean.hat, m),nrow=length(mean.hat)))
  #Get FPC score
  scores <- matrix(NA, nrow=m, ncol=d0)
  for(i in 1:m){
    for (l in 1:d0){
      scores[i,l] = sum(demeaned[i,]*efns0[,l])/n #**
    }
  }
  #get variance structure 
  covb<-array(0,dim = d0)
  for (l in 1:d0){
    for(i in 1:m){
      covb[l] = covb[l] + sum((x[i,]-mean.hat)*efns0[,l])*sum((x[i,]-mean.hat)*efns0[,l])/(n^2)#**
    }
    covb[l] = covb[l]/m
  }
  for(i in 1:d0){
    if (sum(efns0[1:10,i])<0){
      efns0[,i]<-efns0[,i]*-1
    }
  }
  #output lsit includes: 
  #1.dimension, 2.template profile, 3.eigenvalue, 4.eigenfucntions, 5.fpc-scores and 6.variacne
  list(d=d0,mua=mean.hat,evalue = eval0,efns = efns0,scores = scores,covb = covb)
}


