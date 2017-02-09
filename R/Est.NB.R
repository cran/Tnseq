Est.NB <-
function(countsTable, Gene.id, condition, replicates=replicates){
  grp=unique(condition)
  grp=grp[order(grp)] 
  cds <- newCountDataSet(countsTable, condition)
  cds <- estimateSizeFactors( cds )
  if (replicates){
     cds <- estimateDispersions(cds,fitType="local")
  }else cds <- estimateDispersions(cds, method="blind", sharingMode="fit-only", fitType="local")

  res.i <- nbinomTestLess( cds, grp[1], grp[2] )  # nbinomTestLess is for  grp[1]< grp[2]
  
  # p value per insertion (one-sided)
  pBLess.DES=res.i$pval
  
  est=cbind.data.frame(ID=Gene.id, MeanA=res.i$baseMeanA,MeanB=res.i$baseMeanB, pBLess.DES=pBLess.DES)

  return(est=est)
}
