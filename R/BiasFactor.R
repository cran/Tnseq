BiasFactor <-
function(Location, Count, window=10000){
  
  totalLen=diff(range(Location,na.rm=TRUE))
  n.window= ceiling(totalLen/window)
  w.id=cut(Location, n.window,label=FALSE) # each location gets an id
  SummaryStat=aggregate(Count, by = list(w.id), FUN = "sum") ## if there is no count in a window, this window is deleted from the fitting
  
  w=SummaryStat[,1]
  yw=SummaryStat[,2]
  
  #Spline fit (take log count + 1 to deal with zero count)
  loesfit <- loess(log(yw+1) ~ w)  
  
  # make predictions for each locaction
  ypred=predict(loesfit, data.frame(w=w.id)) # impute for each loc even there is zero count in the window
  yfit=exp(ypred-1)
  bias.factor=yfit/exp(mean(log(yfit))) # prod(bias.factor)=1
  
  # make predictions for each window
  ywpred=predict(loesfit, data.frame(w=w))
  ywfit=exp(ywpred-1)
  
  list(bias.factor=bias.factor,yfit=yfit, w=w, yw=yw, ywfit=ywfit)
}
