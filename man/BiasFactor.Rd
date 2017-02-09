\name{BiasFactor}
\alias{BiasFactor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculate the bias factor induced by the genomic replication process for each insertion}
\description{
The whole genome is first divided into N consecutive regions (windows). The number of read counts is summed
within each window and is modelled as a function of the window using the locally weighed scatterplot smoothing (LOESS) function. 
Finally, we calculate a fitted value for each insertion, and then normalize these fitted values to have the product of these values equal to one. 
The normalized values are the bias factors(one for each insertion). }

\usage{
BiasFactor(Location, Count, window=10000)
}
\arguments{
  \item{Location}{a numeric vector specifying insertion locations.}
  \item{Count}{a numeric vector of read counts over locations.}
  \item{window}{the number of regions to divide the whole genome.}
}
  
\value{
 \item{bias.factor}{bias factor for each insertion.}
 \item{yfit}{the fitted count value for each insertion.}
 \item{w}{window id}
 \item{yw}{total read counts over windows.}
 \item{ywfit}{the fitted count value for each window.}
}

 
\author{
Lili Zhao \email{zhaolili@umich.edu}
}


\examples{

\dontrun{
data(serratia)
countData=serratia[,-c(1,2,3)]
location=serratia$Loc


# obtain bias factor for the first sample
foo=BiasFactor(location, countData[,1],  window=10000)  

# obtain the adjusted counts after the bias correction
countData.adj=countData[,1]/foo$bias.factor


}
} 


