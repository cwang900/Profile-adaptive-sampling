library(splines)
library(MASS)
source("MFPCA_s.R")
source("dpca_i.R")
source("Get dpcaAi_delta.R")
source("Getting ARL_delta.R")

#Function to generate simulated data set
#Both Model 1 and Model 2 are generate from below 
Gen_X<-function(m){
  l=seq(1,20,length.out = 20)
  n=length(l)
  Sigma=matrix(data=NA, nrow=n, ncol=n)
  for(i in 1:n){
    for (j in 1:n){
      Sigma[i,j]= 0.8^abs(i-j)
      #Sigma[i,j]= exp(-abs(i-j)^2/20)
    }
  }
  #selected B-spline basis
  vk <-matrix(bs(1:50,df = 6,degree = 2),nrow = 50,ncol = 6)
  X <-array(0,dim = c(m,50,20))
  for (i in 1:m){
    score_k <- mvrnorm(6,mu = rep(0,20),Sigma = Sigma)
    ksi_k<-score_k
    X[i,,] = vk%*%ksi_k
  }
  list(data = X,basis = vk) #basis is also attached 
}
#X_train<-Gen_X(20000)$data
#X_param <-Gen_X(200)$data

#number of sensor
p <-dim(X_param)[3]  
umin <- 1.5
dpca_list <- vector(mode = "list", length = p)
for (i in 1:p){
  xi <- X_param[,,i]
  #store parameter list of each sensor individually into a list
  dpca_list[[i]] <- dpca_est(xi,d0 = 6) 
}
#scaling process
E <-rep(0,length = p)
for(j in 1:p){
  Covb = sum(dpca_list[[j]]$covb)
  E[j] <- 1/sqrt(Covb)
}

#Get the fpca result for rep of 20
fpca_ans<-matrix(rep(0,length = 42),nrow = 42,ncol = 1)
fpca_fillans<-function(Train,obnum){
  Init_d <-ARLi_initial(Train,obnum)
  Initial<-Init_d$Initial
  stat_hist <- df_hist(Train,Initial,obnum)
  #using rule of 1/200 in equation 18
  L <- quantile(stat_hist,199/200)  
  oc_run04 <-rep(0,length = 20)
  oc_run07 <-rep(0,length = 20)
  oc_run10 <-rep(0,length = 20)
  oc_run14 <-rep(0,length = 20)
  oc_run21 <-rep(0,length = 20)
  oc_run28 <-rep(0,length = 20)
  oc_run42 <-rep(0,length = 20)
  #Parital Shift in sensor 3,4,5,8,13
  #20 can be extended to 500, as we did in simulation, or can be set to be arbitrary 
  #large to achieve more converged result 
  for(i in 1:20){ 
    X_ic<-Gen_X(500)$data
    #selection of mean-shift
    for(mu_degree in c(0.4,0.7,1.0,1.4,2.1,2.8,4.2)){
      X_oc <- X_ic #refresh X_oc
      X_oc[,11:40,3] <- X_oc[,11:40,3]+mu_degree*array(1,dim = c(dim(X_ic)[1],30))
      X_oc[,11:40,4] <- X_oc[,11:40,4]+mu_degree*array(1,dim = c(dim(X_ic)[1],30))
      X_oc[,11:40,5] <- X_oc[,11:40,5]+mu_degree*array(1,dim = c(dim(X_ic)[1],30))
      X_oc[,11:40,8] <- X_oc[,11:40,8]+mu_degree*array(1,dim = c(dim(X_ic)[1],30))
      X_oc[,11:40,13] <- X_oc[,11:40,13]+mu_degree*array(1,dim = c(dim(X_ic)[1],30))
      if (mu_degree == 0.4){oc_run04[i] <- Get_ARLi(X_oc,L,obnum)$RL}
      else if(mu_degree == 0.7){oc_run07[i] <- Get_ARLi(X_oc,L,obnum)$RL}
      else if(mu_degree == 1.0){oc_run10[i] <- Get_ARLi(X_oc,L,obnum)$RL}
      else if(mu_degree == 1.4){oc_run14[i] <- Get_ARLi(X_oc,L,obnum)$RL}
      else if(mu_degree == 2.1){oc_run21[i] <- Get_ARLi(X_oc,L,obnum)$RL}
      else if(mu_degree == 2.8){oc_run28[i] <- Get_ARLi(X_oc,L,obnum)$RL}
      else {oc_run42[i] <- Get_ARLi(X_oc,L,obnum)$RL}
    }
  }
  ansarray<-rep(0,14)
  ansarray[1]<-mean(oc_run04,na.rm=TRUE); ansarray[2]<-sd(oc_run04,na.rm = TRUE) 
  ansarray[3]<-mean(oc_run07,na.rm=TRUE) ;ansarray[4]<-sd(oc_run07,na.rm = TRUE) 
  ansarray[5]<-mean(oc_run10,na.rm=TRUE);ansarray[6]<- sd(oc_run10,na.rm = TRUE) 
  ansarray[7]<-mean(oc_run14,na.rm=TRUE);ansarray[8]<-sd(oc_run14,na.rm = TRUE) 
  ansarray[9]<-mean(oc_run21,na.rm=TRUE);ansarray[10]<-sd(oc_run21,na.rm = TRUE) 
  ansarray[11]<-mean(oc_run28,na.rm=TRUE);ansarray[12]<-sd(oc_run28,na.rm = TRUE) 
  ansarray[13]<-mean(oc_run42,na.rm=TRUE);ansarray[14]<-sd(oc_run42,na.rm = TRUE) 
  ansarray
}
#observation capability = 8,4,1
for(r in c(8,4,1)){ 
  if (r == 8){fpca_ans[1:14,1]<-fpca_fillans(X_train,r)}
  else if (r == 4){fpca_ans[15:28,1]<-fpca_fillans(X_train,r)}
  else {fpca_ans[29:42,1]<-fpca_fillans(X_train,r)}
}


#fill table for fpca_del, del = 0.01,0.1,0.5
#20 can be extended to 500, as we did in simulation, or can be set to aribitrary 
#large to more converged result 
del_ans<-matrix(rep(0,length = 126),nrow = 42,ncol =3)
del_fillans<-function(Train,obnum,del){
  Init_d <-ARLi_init_delta(Train,obnum,del)
  Initial<-Init_d$table
  stat_hist <- df_hist_delta(Train,Initial,obnum,del)
  L <- quantile(stat_hist,199/200)  #using rule of 1/200 118
  oc_run04 <-rep(0,length = 20)
  oc_run07 <-rep(0,length = 20)
  oc_run1 <-rep(0,length = 20)
  oc_run14 <-rep(0,length = 20)
  oc_run21 <-rep(0,length = 20)
  oc_run28 <-rep(0,length = 20)
  oc_run42 <-rep(0,length = 20)
  for(i in 1:20){
    X_ic<-Gen_X(500)$data 
    for(mu_degree in c(0.4,0.7,1.0,1.4,2.1,2.8,4.2)){#levels of mean-shift 
      X_oc <- X_ic
      X_oc[,11:40,3] <- X_oc[,11:40,3]+mu_degree*array(1,dim = c(dim(X_ic)[1],30))
      X_oc[,11:40,4] <- X_oc[,11:40,4]+mu_degree*array(1,dim = c(dim(X_ic)[1],30))
      X_oc[,11:40,5] <- X_oc[,11:40,5]+mu_degree*array(1,dim = c(dim(X_ic)[1],30))
      X_oc[,11:40,8] <- X_oc[,11:40,8]+mu_degree*array(1,dim = c(dim(X_ic)[1],30))
      X_oc[,11:40,13] <- X_oc[,11:40,13]+mu_degree*array(1,dim = c(dim(X_ic)[1],30))
      
      if (mu_degree == 0.4){oc_run04[i] <- Get_ARLi_delta(X_oc,L,obnum,del)$RL}
      else if(mu_degree == 0.7){oc_run07[i] <- Get_ARLi_delta(X_oc,L,obnum,del)$RL}
      else if(mu_degree == 1.0){oc_run1[i] <- Get_ARLi_delta(X_oc,L,obnum,del)$RL}
      else if(mu_degree == 1.4){oc_run14[i] <- Get_ARLi_delta(X_oc,L,obnum,del)$RL}
      else if(mu_degree == 2.1){oc_run21[i] <- Get_ARLi_delta(X_oc,L,obnum,del)$RL}
      else if(mu_degree == 2.8){oc_run28[i] <- Get_ARLi_delta(X_oc,L,obnum,del)$RL}
      else {oc_run42[i] <- Get_ARLi_delta(X_oc,L,obnum,del)$RL}
    }
  }
  ansarray<-rep(0,14)
  ansarray[1]<-mean(oc_run04,na.rm=TRUE); ansarray[2]<-sd(oc_run04,na.rm = TRUE) 
  ansarray[3]<-mean(oc_run07,na.rm=TRUE) ;ansarray[4]<-sd(oc_run07,na.rm = TRUE) 
  ansarray[5]<-mean(oc_run1,na.rm=TRUE);ansarray[6]<- sd(oc_run1,na.rm = TRUE) 
  ansarray[7]<-mean(oc_run14,na.rm=TRUE);ansarray[8]<-sd(oc_run14,na.rm = TRUE) 
  ansarray[9]<-mean(oc_run21,na.rm=TRUE);ansarray[10]<-sd(oc_run21,na.rm = TRUE) 
  ansarray[11]<-mean(oc_run28,na.rm=TRUE);ansarray[12]<-sd(oc_run28,na.rm = TRUE) 
  ansarray[13]<-mean(oc_run42,na.rm=TRUE);ansarray[14]<-sd(oc_run42,na.rm = TRUE) 
  ansarray
}

#fill table for fpca_delta, del = 0.01,0.1,0.5
#observation capability is r = 8,4,1
for (del in c(0.01,0.1,0.5)){  #choice of delta
  if (del == 0.01){
    for(r in c(8,4,1)){
      if (r == 8){del_ans[1:14,1]<-del_fillans(X_train,r,del)}
      else if (r == 4){del_ans[15:28,1]<-del_fillans(X_train,r,del)}
      else {del_ans[29:42,1]<-del_fillans(X_train,r,del)}
    }
  }
  else if (del == 0.1){
    for(r in c(8,4,1)){
      if (r == 8){del_ans[1:14,2]<-del_fillans(X_train,r,del)}
      else if (r == 4){del_ans[15:28,2]<-del_fillans(X_train,r,del)}
      else {del_ans[29:42,2]<-del_fillans(X_train,r,del)}
    }
  }
  else { #del == 0.5
    for(r in c(8,4,1)){
      if (r == 8){del_ans[1:14,3]<-del_fillans(X_train,r,del)}
      else if (r == 4){del_ans[15:28,3]<-del_fillans(X_train,r,del)}
      else {del_ans[29:42,3]<-del_fillans(X_train,r,del)}
    }
  }
}
#mfpca parameter estimation 
mf_est <- mfpca_est(X_param,6) #use dimension of 6 eigenvector

Covb_mfpca <- array(0,dim = c(p,p))
#aggregate all the d = 1,2,...d0 for the covariance matrix 
for (i in 1:mf_est$d){
  Covb_mfpca = Covb_mfpca + mf_est$covb[i,,]
}
Corr <- array(0,dim = c(p,p))
E1 <- solve(diag(sqrt(diag(Covb_mfpca)))) 
#Estimated correlation matrix after scaling
Corr <- E1%*%Covb_mfpca%*%E1

mfpca_ans<-matrix(rep(0,length = 42),nrow = 42,ncol = 1)
mfpca_fillans<-function(Train,Corr,obnum){
  Init_d <-ARL_initial(Train,obnum)
  Initial<-Init_d$table
  stat_hist <- mf_hist(Train,Initial,Corr,obnum)
  L <- quantile(stat_hist,199/200)  #using rule of 1/200 118
  oc_run04 <-rep(0,length = 20)
  oc_run07 <-rep(0,length = 20)
  oc_run10 <-rep(0,length = 20)
  oc_run14 <-rep(0,length = 20)
  oc_run21 <-rep(0,length = 20)
  oc_run28 <-rep(0,length = 20)
  oc_run42 <-rep(0,length = 20)
  #Entire Shift
  for(i in 1:20){
    X_ic<-Gen_X(500)$data
    for(mu_degree in c(0.4,0.7,1.0,1.4,2.1,2.8,4.2)){
      X_oc <- X_ic #refresh X_oc
      X_oc[,11:40,3] <- X_oc[,11:40,3]+mu_degree*array(1,dim = c(dim(X_ic)[1],30))
      X_oc[,11:40,4] <- X_oc[,11:40,4]+mu_degree*array(1,dim = c(dim(X_ic)[1],30))
      X_oc[,11:40,5] <- X_oc[,11:40,5]+mu_degree*array(1,dim = c(dim(X_ic)[1],30))
      X_oc[,11:40,8] <- X_oc[,11:40,8]+mu_degree*array(1,dim = c(dim(X_ic)[1],30))
      X_oc[,11:40,13] <- X_oc[,11:40,13]+mu_degree*array(1,dim = c(dim(X_ic)[1],30))
      if (mu_degree == 0.4){oc_run04[i] <- Get_ARL(X_oc,Corr,L,obnum)$RL}
      else if(mu_degree == 0.7){oc_run07[i] <- Get_ARL(X_oc,Corr,L,obnum)$RL}
      else if(mu_degree == 1.0){oc_run10[i] <- Get_ARL(X_oc,Corr,L,obnum)$RL}
      else if(mu_degree == 1.4){oc_run14[i] <- Get_ARL(X_oc,Corr,L,obnum)$RL}
      else if(mu_degree == 2.1){oc_run21[i] <- Get_ARL(X_oc,Corr,L,obnum)$RL}
      else if(mu_degree == 2.8){oc_run28[i] <- Get_ARL(X_oc,Corr,L,obnum)$RL}
      else {oc_run42[i] <- Get_ARL(X_oc,Corr,L,obnum)$RL}
    }
  }
  ansarray<-rep(0,14)
  ansarray[1]<-mean(oc_run04,na.rm=TRUE); ansarray[2]<-sd(oc_run04,na.rm = TRUE) 
  ansarray[3]<-mean(oc_run07,na.rm=TRUE) ;ansarray[4]<-sd(oc_run07,na.rm = TRUE) 
  ansarray[5]<-mean(oc_run10,na.rm=TRUE);ansarray[6]<- sd(oc_run10,na.rm = TRUE) 
  ansarray[7]<-mean(oc_run14,na.rm=TRUE);ansarray[8]<-sd(oc_run14,na.rm = TRUE) 
  ansarray[9]<-mean(oc_run21,na.rm=TRUE);ansarray[10]<-sd(oc_run21,na.rm = TRUE) 
  ansarray[11]<-mean(oc_run28,na.rm=TRUE);ansarray[12]<-sd(oc_run28,na.rm = TRUE) 
  ansarray[13]<-mean(oc_run42,na.rm=TRUE);ansarray[14]<-sd(oc_run42,na.rm = TRUE) 
  ansarray
}
#fill table for fpca
for(r in c(8,4,1)){
  if (r == 8){mfpca_ans[1:14,1]<-mfpca_fillans(X_train,Corr,r)}
  else if (r == 4){mfpca_ans[15:28,1]<-mfpca_fillans(X_train,Corr,r)}
  else {mfpca_ans[29:42,1]<-mfpca_fillans(X_train,Corr,r)}
}

#merge three method 
output = matrix(rep(0,length = 210),ncol = 5)
output[,1] = mfpca_ans[,1]
output[,2] = fpca_ans[,1]
output[,3:5] = del_ans[,1:3]

#show just the ARL, but not SDRL
output_nosd<-output[c(1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39,41),]
out_table<-as.data.frame(output_nosd)
colnames(out_table)<-c('proposed','fpca','del05','del01','del001')
out_table
