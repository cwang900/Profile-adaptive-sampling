#FPCA + N(0,1) approach CUSUM procedure

#set up initial table
ARLi_initial<-function(testset,obnum){
r<-obnum
cusum_score<-rep(0,length=p)
for (i in 1:r){
  ob_init<-dpca_score(testset[1,,i],dpca_list[[i]]$evalue,dpca_list[[i]]$efns,dpca_list[[i]]$mua) #d*1
  cusum_score[i] <- sum(ob_init)*E[i] #stanardlized FPC score 
}

for(i in (r+1):p){
  ran_init <- rnorm(1,0,1) ## i.i.d N(0,1) 
  cusum_score[i] <- ran_init  
}
#Two sided CUSUM W_i = max(W_pos,W_neg)
Initial <- data.frame(sensor_id = seq(1,20,length.out = 20),cusum_score)
for(i in 1:length(cusum_score)){
  Initial[i,"W_positive"]<-max(umin*cusum_score[i]-(umin^2)/2,0)
  Initial[i,"W_negative"]<-max(-umin*cusum_score[i]-(umin^2)/2,0)
  Initial[i,"W_i"]<-max(Initial[i,"W_positive"],Initial[i,"W_negative"])
}
#reorder local statistics
Initial["rank"] = rank(-Initial["W_i"],ties.method = "random")
ob_sensor <- subset(Initial,rank<= r)$sensor_id
list(Initial = Initial,observed = ob_sensor)
}

#Get Running Length of 1 repetition 
#L is the control limits, obnum is the observation number 
Get_ARLi<-function(testset,L,obnum){
  r<-obnum
  procedure <-vector(mode ="list",length = 500)
  ob_hist<-array(0,dim = c(500,r))
  L_hist <- rep(0,dim(testset)[1]) #depend on train size
  init_list <- ARLi_initial(testset,obnum)
  ob_hist[1,] <- ob_hist[1,]+ init_list$observed
  procedure[[1]]<-init_list$Initial
  s <- as.matrix(init_list$Initial["W_i"])
  #Global index 
  L_hist[1] <- sqrt(t(s)%*%solve(diag(length(s)))%*%s)  
  if(L_hist[1]>=L){
    return (list(RL = 1,observed = ob_hist,procedure = procedure,L_trace = L_hist))
  }
  current <- init_list$Initial
  for(sample in 2:dim(testset)[1]){
    next_step <- cusum_df(testset,current,sample,obnum) #return the latest table
    next_table<- next_step$table
    next_ob<-next_step$observed 
    s <- as.matrix(next_table["W_i"])
    #Global Index 
    L_hist[sample] <- sqrt(t(s)%*%solve(diag(length(s)))%*%s) ##think about it
    procedure[[sample]]<-next_table
    ob_hist[sample,] <- ob_hist[sample,] + next_ob
    #Check if exceed the control limits 
    if(L_hist[sample]>=L){
      return(list(RL = sample,ob_hist = ob_hist,procedure = procedure,L_trace = L_hist))
    }  
    current <- next_table
  }
  return(list(RL = 500,ob_hist = ob_hist,procedure = procedure,L_hist = L_hist))
}

#1.get histgram for control limit L to apply 1/200 rule
df_hist<-function(Entire,Initial,obnum){
  L_hist <- rep(0,dim(Entire)[1]) #depend on train size
  s <- as.matrix(Initial["W_i"])
  L_hist[1] <- sqrt(t(s)%*%solve(diag(length(s)))%*%s) #think about it 
  current <- Initial
  for(sample in 2:dim(Entire)[1]){
    next_step <- cusum_df(Entire,current,sample,obnum) 
    next_table<-next_step$table #return the latest table
    s <- as.matrix(next_table["W_i"])
    L_hist[sample] <- sqrt(t(s)%*%solve(diag(length(s)))%*%s)
    current <- next_table
  }
  #return hist of control limits 
  L_hist
}

#CUSUM procedure and sensor redistribution for FPCA
cusum_df <- function(Entire,input,example,obnum){
  r<-obnum
  #Might need to think of truncate the matrix
  table <- input
  for (i in 1:p){
    if(table[i,]$rank<=r){
      #Compensation for the observed case
      scale = E[i]
      ob_comp <- dpca_score(Entire[example,,i],dpca_list[[i]]$evalue,dpca_list[[i]]$efns,dpca_list[[i]]$mua) #d*1
      table[i,"cusum_score"] <- sum(ob_comp)*scale #N(0,1) variable
    }
    else{
      #compensation for the unobserved case, which comes from individual N(0,1)
      ran_comp <- rnorm(1,0,1) #N(0,1) normal variable 
      table[i,"cusum_score"] <-ran_comp
    }
    #Local CUSUM
    table[i,"W_positive"]<-max(table[i,"W_positive"]+umin*table[i,"cusum_score"]-(umin^2)/2,0)
    table[i,"W_negative"]<-max(table[i,"W_negative"]-umin*table[i,"cusum_score"]-(umin^2)/2,0)
    table[i,"W_i"]<-max(table[i,"W_positive"],table[i,"W_negative"])
  }
  table$rank = rank(-1*table$W_i,ties.method = "random")
  ob_sensor <- subset(table,rank<=r)$sensor_id
  list(table = table,observed = ob_sensor)
}



#FPCA + Del = 0.01,0.1,0.5 approach under CUSUM framework

ARLi_init_delta<-function(testset,obnum,delta){
  r<-obnum
  cusum_score<-rep(0,length=p)
  for (i in 1:r){
    ob_init<-dpca_score(testset[1,,i],dpca_list[[i]]$evalue,dpca_list[[i]]$efns,dpca_list[[i]]$mua) 
    cusum_score[i] <- sum(ob_init)*E[i]
  }
  for(i in (r+1):p){
    #unobserved compensation become delta = 0.01,0.1,0.5
    cusum_score[i] <- delta 
  }
  #update two sided local CUSUM
  table <- data.frame(sensor_id = seq(1,20,length.out = 20),cusum_score)
  for(i in 1:length(cusum_score)){
    table[i,"W_positive"]<-max(umin*cusum_score[i]-(umin^2)/2,0)
    table[i,"W_negative"]<-max(-umin*cusum_score[i]-(umin^2)/2,0)
    table[i,"W_i"]<-max(table[i,"W_positive"],table[i,"W_negative"])
  }
  table["rank"] = rank(-table["W_i"],ties.method = "random")
  ob_sensor <- subset(table,rank<= r)$sensor_id
  list(table = table ,observed = ob_sensor)
}


Get_ARLi_delta<-function(testset,L,obnum,delta){
  procedure <-vector(mode ="list",length = 500)
  ob_hist<-array(0,dim = c(500,r))
  L_hist <- rep(0,dim(testset)[1]) #depend on train size
  init_list <- ARLi_init_delta(testset,obnum,delta)
  ob_hist[1,] <- ob_hist[1,]+ init_list$observed
  procedure[[1]]<-init_list$table
  s <- as.matrix(init_list$table["W_i"])
  L_hist[1] <- sqrt(t(s)%*%solve(diag(length(s)))%*%s)  
  if(L_hist[1]>=L){
    return (list(RL = 1,observed = ob_hist,procedure = procedure,L_trace = L_hist))
  }
  current <- init_list$table
  for(sample in 2:dim(testset)[1]){
    next_step <- cusum_df_delta(testset,current,sample,obnum,delta) 
    next_table<- next_step$table  
    next_ob<-next_step$observed 
    s <- as.matrix(next_table["W_i"])
    L_hist[sample] <- sqrt(t(s)%*%solve(diag(length(s)))%*%s)
    procedure[[sample]]<-next_table
    ob_hist[sample,] <- ob_hist[sample,] + next_ob
    if(L_hist[sample]>=L){
      return(list(RL = sample,ob_hist = ob_hist,procedure = procedure,L_trace = L_hist))
    }  
    current <- next_table
  }
  return(list(RL = 500,ob_hist = ob_hist,procedure = procedure,L_hist = L_hist))
}

df_hist_delta<-function(Entire,Initial,obnum,delta){
  L_hist <- rep(0,dim(Entire)[1])
  s <- as.matrix(Initial["W_i"])
  L_hist[1] <- sqrt(t(s)%*%solve(diag(length(s)))%*%s) 
  current <- Initial
  for(sample in 2:dim(Entire)[1]){
    next_step <- cusum_df_delta(Entire,current,sample,obnum,delta) 
    next_table<-next_step$table 
    s <- as.matrix(next_table["W_i"])
    L_hist[sample] <- sqrt(t(s)%*%solve(diag(length(s)))%*%s)
    current <- next_table
  }
  L_hist
}

cusum_df_delta <- function(Entire,input,example,obnum,delta){
  r<-obnum
  table <- input
  for (i in 1:p){
    if(table[i,]$rank<=r){
      scale = E[i]
      ob_comp <- dpca_score(Entire[example,,i],dpca_list[[i]]$evalue,dpca_list[[i]]$efns,dpca_list[[i]]$mua) #d*1
      table[i,"cusum_score"] <- sum(ob_comp)*scale #N(0,1) variable
    }
    else{
      #unobserved compensation become delta = 0.01,0.1,0.5 from FPCA
      table[i,"cusum_score"] <-delta
    }
    table[i,"W_positive"]<-max(table[i,"W_positive"]+umin*table[i,"cusum_score"]-(umin^2)/2,0)
    table[i,"W_negative"]<-max(table[i,"W_negative"]-umin*table[i,"cusum_score"]-(umin^2)/2,0)
    table[i,"W_i"]<-max(table[i,"W_positive"],table[i,"W_negative"])
  }
  table$rank = rank(-1*table$W_i,ties.method = "random")
  ob_sensor <- subset(table,rank<=r)$sensor_id
  list(table = table,observed = ob_sensor)
}


