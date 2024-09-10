#X1 will be a single profile sample readings, so dim(X1) = 1*n*p 
#eval,efns,mua are from estimated mfpca model 
m_score<-function(X1,eval,efns,mua){
  demeaned = X1-mua[1,,] #n*p
  n <- dim(efns)[1]
  p <- dim(mua)[3]
  #scores
  d = length(eval)
  scores <-array(0,dim = c(d,p))
  #see E.q 10
  for(l in 1:d){
    scores[l,] = colSums(demeaned*efns[,l])/n 
  }
  scores #d*p matrix
}

#Input X is of m*n*p matrix (number of samples, number of time points, number of sensors)
##d0 is the number of decomposition dimension. in E.q 5
mfpca_est<-function(X,d0){
  m = dim(X)[1]
  n = dim(X)[2]
  p = dim(X)[3]
  mua <- array(0,dim = c(1,n,p))
  # Equation 6 in the paper 
  for(i in 1:p){
    for(j in 1:n){
      mua[1,j,i] = mean(X[,j,i]) 
    }
  }
  #Convarinace function in E.q 7 in the paper 
  cova = array(0,dim = c(n,n))
  for(u in 1:m){
    sum_p = array(0,dim = c(n,n))
    for (i in 1:p){
      sum_p = sum_p + (X[u,,i]-mua[1,,i])%*%t(X[u,,i]-mua[1,,i]) 
    }
    cova = cova + sum_p
  }
  cova = cova/m
  #Eigen Decompostion is in E.q 8  
  eio<-eigen(cova)
  eval<-eio$values/n
  evec<-eio$vectors*sqrt(n)
  #d0 = 6
  d = d0
  eval0 <- eval[1:d]
  efns0 <- evec[,1:d]
  demeaned <- array(0,dim = c(m,n,p))
  for (i in 1:m){
    demeaned[i,,] = X[i,,]-mua[1,,]
  }
  #scores that is calculated by E.q 10
  scores <-array(0,dim = c(m,d,p))
  for(i in 1:m){
    for(l in 1:d){
      scores[i,l,] = colSums(demeaned[i,,]*efns0[,l])/n
    }
  }
  #Estimation covariance matrix of MFPCA score Ksi, introduced in E.q 9
  covb<-array(0,dim=c(d,p,p))
  for(l in 1:d){
    for (j in 1:p){
      for (h in 1:p){
        for(i in 1:m){
          a = (X[i,,j]-mua[1,,j])*efns0[,l]/n
          b = (X[i,,h]-mua[1,,h])*efns0[,l]/n
          sum_n = sum(a)*sum(b)
          covb[l,j,h] = covb[l,j,h]+sum_n
        }
      }
    } 
  }
  covb = covb/m
  #Adjust sign of efns 
  for(i in 1:d){
    if (sum(efns0[1:10,i])<0){ #using the second time-stamp 
      efns0[,i]<-efns0[,i]*-1
    }
  }
  #output lsit includes: 
  #1.dimension, 2.template profile, 3.eigenvalue, 4.eigenfucntions, 5.mfpc-scores and 6.covariacne matrix 
  list(d=d0,mua=mua,evalue = eval0,efns = efns0,scores = scores,covb = covb)
}
