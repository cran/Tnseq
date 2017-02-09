CDfun.exact <-
function (dlist=dlist, limit=limit){
  
  # logFC is the median logFC
  logFC.cd=sapply(dlist,function(l) 
    uniroot(function(x)  sum(qnorm((pt((x-l$bhat)/l$se,l$df)))), 
            lower = -limit, upper = limit, trace=TRUE, tol = 1e-10)$root
  )
  
  pvalBLess.cd=sapply(dlist, function(l)  1-pnorm(1/sqrt(nrow(l))*sum(qnorm(pt(-l$bhat/l$se,l$df)))))
  pvalALess.cd= 1- pvalBLess.cd
  
  pval.cd= 2*pmin(pvalBLess.cd,pvalALess.cd)

  res.cd=cbind.data.frame(ID=names(pvalBLess.cd), FC=2^(logFC.cd),log2FC=logFC.cd,pvalue=pval.cd)
  
  return(res.cd)  
}
