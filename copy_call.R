copy_call = function(rid_obj,knn=0,w_size=5,prob_thresh=.025,mult=c(0,.5,1,1.5,2,3)){
  library(stringr)
  library(MASS)
  library(fitdistrplus)

  knn = knn+1

  #load target window
  wnd = rid_obj$window
  Yids = wnd[,1] == "chrY" | wnd[,1] == "chrX"

  #Load model residuals
  fv = rid_obj$expected
  nv = rid_obj$obs_norm
  fv = fv[!Yids,]
  nv = nv[!Yids,]
  wnd = wnd[!Yids,]

  #Prune low coverage windows at end of chromosomes
  for(chr in unique(wnd[,1])){
    cids = which(wnd[,1] == chr)
    lid = cids[length(cids)]
    ratio = (wnd[lid,3]-wnd[lid,2])/(wnd[lid-1,3]-wnd[lid-1,2])
    if(ratio < .75){
      fv[lid,] = nv[lid,]
    }
  }

  print("Data ready")

  #Identify k-nearest neighbors for bulk profiling
  if(knn ==1){
    kids = matrix(1:ncol(nv),nrow=ncol(nv),ncol=knn)
  } else{
    res = log((nv+.1)/(fv+.1))
    res[is.na(res)] = 0
    D = as.matrix(dist(t(res),method="euclidean"))
    kids = matrix(0,nrow=ncol(nv),ncol=knn)
    for(i in 1:nrow(kids)){
      kids[i,] = order(D[i,],decreasing=FALSE)[1:knn]
    }
  }

  print("KNN computed") 

  tot_size = knn*(2*w_size+1)

  #Calulate model copy numbers
  emp_res = 1000
  qs = matrix(0,nrow=emp_res+1,ncol=tot_size)
  for(i in 1:tot_size){
    su = rowMeans(matrix(runif(100000*i),ncol=i))
    qs[,i] = quantile(su,c(0,(1:emp_res)/emp_res))
  }
  rv = matrix(1,nrow=nrow(fv),ncol=ncol(fv))
  colnames(rv) = colnames(fv)
  rownames(rv) = wnd[,4]
  pv = .5*rv
  breaks = mult
  mult[mult==0] = .1
  
  print("Empirical distribution ready") 
  
  for(i in 1:nrow(nv)){
    chr = wnd[i,1]
    ids = (i-w_size):(i+w_size)
    ids = ids[ids > 0 & ids <= nrow(nv)]
    ids = ids[wnd[ids,1]==wnd[i,1]]
    for(j in 1:ncol(nv)){
      if(nv[i,j]== -1){next}
      nv_j = as.vector(nv[ids,kids[j,]])
      fv_j = as.vector(fv[ids,kids[j,]])
      fv_j = fv_j[nv_j > -1]
      nv_j = nv_j[nv_j > -1]
      mv = mean(ppois(nv_j,fv_j))
      
      num = length(nv_j)
      pv[i,j] = (which.min(abs(qs[,num]-mv))-1)/emp_res

      LL = rep(0,length(mult))
      for(k in 1:length(mult)){
        LL[k] = sum(dpois(nv_j,mult[k]*fv_j,log=TRUE))        
      }
      rv[i,j] = breaks[which.max(LL)]
      if(is.infinite(max(LL))){
        if(mv > .5){rv[i,j] = max(breaks)}
        if(mv < .5){rv[i,j] = min(breaks)}
      }
    }
  }

  print("Probabilities calculated")
  
  pvl = c(pv)
  pvh = c(pv)

  ids = unlist(nv) != -1

  pvl[ids] = p.adjust(pvl[ids],method="fdr")
  pvh[ids] = p.adjust((1-pvh)[ids],method="fdr")

  print("FDR correction complete")

  pvl = matrix(pvl,nrow=nrow(pv),ncol=ncol(pv))
  pvh = matrix(pvh,nrow=nrow(pv),ncol=ncol(pv))

  rv[pvl > prob_thresh & pvh > prob_thresh] = 1 

  #Add copy calls to object
  rid_obj$CNV = rv
  rid_obj$prob = pv
  return(rid_obj)
}
