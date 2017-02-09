CombinePvals <-
function(nb.est=nb.est){
  p=nb.est$pBLess.DES
  p.other= 1- p
  
  dLBLess=split(p, nb.est$ID)  
  dLALess=split(p.other,nb.est$ID)
  
  pBLess=sapply(dLBLess,function(s) stouffer(pmax(1.0e-8,pmin(s,0.99999))))
  pALess=sapply(dLALess,function(s) stouffer(pmax(1.0e-8,pmin(s,0.99999))))


  pval.nb=2*pmin(pBLess,pALess)  # two-sided p
  
  res.NB=cbind.data.frame(ID=names(pBLess),p.nb=pval.nb)
  return(res.NB)
}
