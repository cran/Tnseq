CombineNBEst <-
function(nb.est=nb.est){
  
  
  res.nb=CombinePvals(nb.est=nb.est)  
    
  # obtain averaged total counts per gene
  TotalCount=aggregate(as.matrix(nb.est[,c(2,3)]) ~ ID, data = nb.est, sum) #nb.est[,c(2,3)is MeanA and MeanB
  
  res=merge(TotalCount,res.nb, by="ID",all=TRUE)
  
  res$FC=res$MeanB/res$MeanA
  res$log2FC=log2(res$FC)
  resc=res[,c("ID","FC","log2FC","p.nb")]
  colnames(resc)[4]="pvalue"
  return(resc)
}
