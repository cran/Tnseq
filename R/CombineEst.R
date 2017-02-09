CombineEst <-
function(limma.est=limma.est, nb.est=nb.est, weights="equal", p.nb=FALSE, limit=limit){
  
  if(weights=="hc"){  
    
    #weigts based on input counts only      
    m1=limma.est$Mean1
    hc=Ckmeans.1d.dp(m1, k=c(1,2)) #A fast dynamic programming algorithm for optimal univariate k-means clustering.
    small.counts=m1[hc$cluster==1]
    min=min(small.counts)
    max=max(small.counts)
    print(paste("Weight<1 if Mean1(log)<", round(max,2),sep=""))
    slope= 3/(max-min)  # 0 = e10-3
    intercept=-max*slope
    ww=exp(intercept+slope*m1)
    limma.est$weights=ifelse(ww>1,1,ww)
  
    dlist=split(limma.est, limma.est$ID)   
    res.cd=CDfun.weights(dlist=dlist,limit=limit)
    
  }else if(weights=="equal"){
    dlist=split(limma.est, limma.est$ID)  
    res.cd=CDfun.exact(dlist=dlist, limit=limit)
    
  }else{    
    if (length(weights) != dim(limma.est)[1]) stop("length of weights is not equal to the number of insertions")
    limma.est$weights=weights
    dlist=split(limma.est, limma.est$ID)   
    res.cd=CDfun.weights(dlist=dlist,limit=limit)
  }
  
  if(p.nb==TRUE){
    res.nb=CombinePvals(nb.est=nb.est)  
    res=merge(res.cd, res.nb,  by=c("ID"), all=TRUE)     
  }else res=res.cd
  
  # obtain averaged total counts per gene
  TotalCount=aggregate(as.matrix(nb.est[,c(2,3)]) ~ ID, data = nb.est, sum) #nb.est[,c(2,3)is MeanA and MeanB
  
  res.all=merge(TotalCount,res,by="ID",all=TRUE)
  
  list(res.all=res.all, est=limma.est)
}
