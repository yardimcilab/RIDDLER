riddler_fit = function(infile,feat_file,peak_file,cov_min=10000){
  library(robustbase)
  library(MASS)
  library(stringr)
  library(fitdistrplus)

  M = read.table(infile,header=TRUE,sep=" ")
  R = M[,-(1:4)]
  wnd = M[,(1:4)]
  coverage = colSums(R)
  cids = coverage >= cov_min
  R = R[,cids]
  M = cbind(M[,1:4],R)
  barcodes = colnames(R)

  X = read.csv(feat_file)
  P = read.table(peak_file,sep=" ")
  X = cbind(X,peaks=P[,5])
  #gc,map,peaks
  X = X[,c("gc","map","peaks")]
  for(i in 1:ncol(X)){
    X[,i] = X[,i] / quantile(X[,i],.75)
  } 
 
  #Calculate window probabilities and coverage thresholds
  Zp = matrix(0,nrow=nrow(M),ncol=ncol(R))
  prob = rowSums(R)
  prob = prob/sum(prob)
  for(i in 1:ncol(Zp)){
    cov = sum(R[,i])
    for(j in 1:nrow(Zp)){
      Zp[j,i] = (1-prob[j])^cov
    }
  }
  Zp = matrix(p.adjust(as.vector(Zp)),nrow=nrow(Zp))
  prob_thresh = .05
    
  n = ncol(R)
  fv = matrix(-1,nrow=nrow(M),ncol=n)
  colnames(fv) = barcodes
  nv = fv
  mask = 0*nv
   
  #quantile normalization
  q_list = c(0,1:100)/100
  mu = rep(0,n)
  size = rep(0,n)
  for(i in 1:n){
    Y = R[,i]
    cov = sum(Y)
    ids = Y > 0 | Zp[,i] <= prob_thresh    
    Y = Y[ids]
      
    scale = mean(Y,trim=.15)
    if(scale==0){ scale = mean(Y,trim=.15) }      
    Y = round(100 * Y / scale)

    nv[ids,i] = Y
    mask[ids,i] = 1
    if(i %% 1000 == 0){
      print(i)
    }
  }

  #Parse features and fit model
  Yt = as.vector(nv[mask==1])
  feat_num = ncol(X)
  Xt = do.call("rbind", rep(list(X), n))
  nids = as.vector(mask==1)
  Xt = Xt[nids,]
  
  model = glmrob(Yt ~ .,data=Xt,family="poisson",control= glmrobMqle.control(tcc=1.345),start=c(1,rep(0.1,ncol(Xt))))
  coef = model$coefficients
  fvals = exp(predict(model,newdata=X))
  for(i in 1:n){
    fv[,i] = fvals
  }
  fv[mask!=1] = -1
  #create and return object
  rid_obj = list(obs_norm=nv,expected=fv,model=model,window=wnd)

  return(rid_obj)
}
