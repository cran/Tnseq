stouffer <-
function(x){
  x <- x[!is.na(x)]
  pnorm(sum(qnorm(x)) / sqrt(length(x)))
}
