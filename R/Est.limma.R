Est.limma <-
function(countsTable, Gene.id, condition, norm=TRUE, bayes=bayes){
  
  design = model.matrix(~factor(condition))  # model effect of B over A ==> beta (e.g logFC) <0 is the alternative
  dge <- DGEList(counts=countsTable)
  if(norm){
    dge <- calcNormFactors(dge)
  }else  dge <- calcNormFactors(dge,method="none")  # no normalization
  v <- voom(dge,design,plot=FALSE)

  # numeric matrix of normalized expression values on the log2 scale
  normDat=v$E
  grp=unique(condition)
  n=table(condition)
  a=matrix(normDat[,condition==grp[1]],ncol=n[1])
  b=matrix(normDat[,condition==grp[2]],ncol=n[2])
  
  # mean over replicates for each insertion
  Mean1=rowMeans(a); Mean2=rowMeans(b)
  
  # fit limma model
  fitlimma=lmFit(v)  
  fitbayes = eBayes(fitlimma)
  if(bayes){
    # Empirical bayes
    est=cbind.data.frame(ID=Gene.id, Mean1=Mean1,Mean2=Mean2,df=fitbayes$df.total,bhat=fitbayes$coefficients[,2], p=fitbayes$p.value[,2], se=(fitbayes$stdev.unscaled*sqrt(fitbayes$s2.post))[,2])
      
  }else{
    # Not Empirical bayes, only voom
    est=cbind.data.frame(ID=Gene.id, Mean1=Mean1,Mean2=Mean2,df=fitlimma$df.residual,bhat=fitlimma$coefficients[,2],  se=(fitlimma$stdev.unscaled*fitlimma$sigma)[,2])
      
  }
 
  return(est=est)
}
