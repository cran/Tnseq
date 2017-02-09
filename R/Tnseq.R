Tnseq<-function (countData, geneID, location, pool, condition, weights = "equal", 
          bayes = FALSE, p.nb = FALSE, norm = TRUE, cut=0) 
{
  
  if (!all.equal(nrow(countData), length(geneID), length(location))) 
    stop("length of geneID/location should be equal to the number\n  of rows in countData")
  if (!all.equal(ncol(countData), length(condition), length(pool))) 
    stop("length of condition/pool should be equal to the number\n  of columns in countData")  
  if (!is.character(geneID)) 
    stop("geneID must be a character string")

  grp = unique(condition)
  grp = grp[order(grp)]
  pid = unique(pool)
  npool = length(pid)
  # check replicates per pool (replicates = FALSE if no replicate in both condition in any pool)
  rep.pool=colSums(table(condition, pool))
  replicates = all(rep.pool>2)

  # calcuate the number of (uniuqe) insertions for each gene under each group
  dd.NumInser=NULL
  for (k in 1:length(grp)){
  id.tmp=loc.tmp=NULL
   for (i in 1:npool) {
    sel=(pool == pid[i] & condition == grp[k])
    CountData.grp =countData[, sel] 
    if(sum(sel)==1){
      nonzero = CountData.grp > 0
    } else nonzero = rowSums(CountData.grp) > 0 
    
    id.tmp=c(id.tmp,as.character(geneID[nonzero]))
    loc.tmp=c(loc.tmp,location[nonzero])
  }
  dd2=cbind.data.frame(id=id.tmp,loc=loc.tmp)
  dd2list=split(dd2,dd2$id)
  NumInser=t(sapply(dd2list, function(s) c(NumInser=length(s$loc),Unique.NumInser=length(unique(s$loc)))))
  dd.NumInser[[k]]=cbind.data.frame(ID=rownames(NumInser),NumInser)
  }

  # prepare data for the analysis
  CountData = GeneID = Location = Condition = NULL
  for (i in 1:npool) {
    countsTable = countData[, pool == pid[i]]
    nonzero = rowSums(countsTable) > cut
    CountData[[i]] = countsTable[nonzero, ]
    GeneID[[i]] = geneID[nonzero]
    Location[[i]] = location[nonzero]
    Condition[[i]] = condition[pool == pid[i]]
  }
  
  limma.est.pool = nb.est.pool = NULL
  for (l in 1:length(Condition)) {
    gene.id = GeneID[[l]]
    count.matrix = CountData[[l]]
    condition = Condition[[l]]
    nb.est.pool[[l]] = Est.NB(count.matrix, gene.id, condition, 
                              replicates = replicates)
    if (replicates) {
      limma.est.pool[[l]] = Est.limma(count.matrix, gene.id, condition, bayes = bayes, norm=TRUE)
    }
    print(paste("Finish Estimation for Pool ", l, sep = ""))
  }
  nb.est = do.call("rbind", nb.est.pool)
  if (replicates) {
    limma.est = do.call("rbind", limma.est.pool)
    results = CombineEst(limma.est = limma.est, nb.est = nb.est, 
                         weights = weights, p.nb = p.nb, limit = 100)
    res = results$res.all
  }
  else {
    res = CombineNBEst(nb.est = nb.est)
  }
  res = merge(dd.NumInser[[1]], res, by = "ID", all.y = TRUE)  # only include insertion data from the first group
  
  names(res) <- gsub("^MeanA$", paste("Mean", grp[1], sep = "."), 
                     names(res))
  names(res) <- gsub("^MeanB$", paste("Mean", grp[2], sep = "."), 
                     names(res))
  names(res) <- gsub("^log2FC$", paste("log2FC", "(", grp[2], 
                                       "/", grp[1], ")", sep = ""), names(res))
  if (replicates) {
    list(resTable = res, est.insertion = results$est)
  }
  else list(resTable = res)
}



