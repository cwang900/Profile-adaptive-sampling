#mainly for proposed method
ARL_initial<-function(testset,obnum){ 
  r<-obnum
  #set up initial table and setting up initial compensation
  cusum_score<-rep(0,length=p) 
  #each score are sum(ksi1+ksi2+...ksi6)
  data_samp = testset[1,,]
  ob_init <- m_score(data_samp,mf_est$evalue,mf_est$efns,mf_est$mua) #d*p matrix score table
  for(i in 1:r){
    scale <- E1[i,i] #in E.q 13
    cusum_score[i] <- sum(ob_init[,i])*scale #N(0,1) variable
  }
  ran_init <- mvrnorm (mu=rep (0, times = p),Sigma = Corr)
  for(i in (r+1):p){
    cusum_score[i] <- ran_init[i] #N(0,1) variable
  }
  #introduced in E.q 14
  table <- data.frame(sensor_id = seq(1,20,length.out = 20),cusum_score)
  for(i in 1:length(cusum_score)){
    table[i,"W_positive"]<-max(umin*cusum_score[i]-(umin^2)/2,0)
    table[i,"W_negative"]<-max(-umin*cusum_score[i]-(umin^2)/2,0)
    table[i,"W_i"]<-max(table[i,"W_positive"],table[i,"W_negative"])
  }
  #order local statistics in E.q 15
  table["rank"] = rank(-table["W_i"],ties.method = "random")
  ob_sensor <- subset(table,rank<= r)$sensor_id
  list(table = table,observed = ob_sensor)
}
Get_ARL <- function(testset,Corr,L,obnum){
  procedure <-vector(mode ="list",length = 500)
  ob_hist<-array(0,dim = c(500,r))
  L_hist <- rep(0,dim(testset)[1]) #depend on train size
  Init <- ARL_initial(testset,obnum)
  ob_hist[1,] <- ob_hist[1,]+ Init$observed
  procedure[[1]]<-Init$table
  s <- as.matrix(Init$table["W_i"]) #put a note
  L_hist[1] <- sqrt(t(s)%*%solve(Corr)%*%s)
  if(L_hist[1]>=L){
    return (list(RL = 1,observed = ob_hist,procedure = procedure,L_trace = L_hist))
  }
  current <- Init$table
  for(sample in 2:dim(testset)[1]){
    next_step <- cusum_stage(testset,current,sample,obnum)
    next_table<- next_step$table
    next_ob<-next_step$observed 
    s <- as.matrix(next_table["W_i"])
    #constuct of global statistics in E.q 16
    L_hist[sample] <- sqrt(t(s)%*%solve(Corr)%*%s)
    procedure[[sample]]<-next_table
    ob_hist[sample,] <- ob_hist[sample,] + next_ob
    #Check if G(i) exceed control limits L in E.q 17
    if(L_hist[sample]>=L){
      return(list(RL = sample,ob_hist = ob_hist,procedure = procedure,L_hist = L_hist))
    }
    current <- next_table
  }
  return(list(RL = 500,ob_hist = ob_hist,procedure = procedure,L_hist = L_hist))
}

#get histgram for mfpca method L 
mf_hist<-function(Entire,Initial,Corr,obnum){
  L_hist <- rep(0,dim(Entire)[1]) #depend on train size
  s <- as.matrix(Initial["W_i"])
  #in E.q 16
  L_hist[1] <- sqrt(t(s)%*%solve(Corr)%*%s)
  current <- Initial
  for(sample in 2:dim(Entire)[1]){
    next_step <- cusum_stage(Entire,current,sample,obnum)$table #return the latest table
    s <- as.matrix(next_step["W_i"])
    #in E.q 16
    L_hist[sample] <- sqrt(t(s)%*%solve(Corr)%*%s)
    current <- next_step
  }
  L_hist
}
#mainly for proposed method 
cusum_stage <- function(Entire,input,example,obnum){
  table <- input
  r<-obnum
  ob_comp <- m_score(Entire[example,,],mf_est$evalue,mf_est$efns,mf_est$mua) #return d*p
  #random compensation 
  ran_comp = mvrnorm (mu=rep (0, times = p),Sigma = Corr)
  for (i in 1:p){
    if(table[i,]$rank<=r){
      #observed sensor and aggregate all the d = 1,2,...d0 for the mfpc scores
      scale = E1[i,i] #scaling in E.q 13
      table[i,"cusum_score"] = sum(ob_comp[,i])*scale #marginal N(0,1) normal variable 
    }
    else{
      #unobserved sensor sample from multivariate normal
      table[i,"cusum_score"] = ran_comp[i] #marginal N(0,1) normal variable 
    }
    #CUSUM procedure in E.q 14
    table[i,"W_positive"]<-max(table[i,"W_positive"]+umin*table[i,"cusum_score"]-(umin^2)/2,0)
    table[i,"W_negative"]<-max(table[i,"W_negative"]-umin*table[i,"cusum_score"]-(umin^2)/2,0)
    table[i,"W_i"]<-max(table[i,"W_positive"],table[i,"W_negative"])
  }
  #order local statistics in E.q 15
  table$rank <- rank(-1*table$W_i,ties.method = "random")
  ob_sensor <- subset(table,rank<=r)$sensor_id
  list(table = table,observed = ob_sensor)
}