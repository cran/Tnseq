\name{TnseqDiff}
\alias{TnseqDiff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Identify conditionally essential genes using high-throughput sequencing data from transposon mutant libraries
}
\description{
It utilizes two steps to estimate the essentiality for each gene in the genome.
First, it construcsts a confidence distribution (CD) function containing evidence of conditional essentiality for each insertion by comparing read counts
of that insertion between conditions. Second, it combines insertion-level CD functions
to infer the essentiality for the corresponding gene. (Zhao et al., 2017).
If no replicate in both conditions in any pool, the insertion-level p-values are calculated from the DESeq (using method=\code{"blind"} and sharingMode=\code{"fit-only"}), 
and then combined using the stouffer method.
}
\usage{
TnseqDiff(countData, geneID, location, pool, condition, weights="equal", bayes=FALSE,
        p.nb=FALSE, norm=TRUE, cut=0)
}

\arguments{
  \item{countData}{a read count matrix. It has n rows (each row corresponding to a unique insertion site) and m columns (each column corresponding to a sample).}
  \item{geneID}{a character string of gene names for the insertion sites in countData.}
  \item{location}{a numeric vector specifying insertion locations.}
  \item{pool}{a numeric vector specifying the pool id that each sample in countData belongs to.}
  \item{condition}{a character string specifying the condition that each sample in countData belongs to.}
  \item{weights}{a character string specifying weights for the insertion sites. These weights are used to weight the insertion-level evidence for the combination. 
    It must be \code{"equal"} (default) or \code{"hc"}; \code{"equal"} assumes equal weights for all insertions, and \code{"hc"} generates small weights for low insertion counts in the input condition. These weights are determined using the algorithm described in (Wang and Song, 2001).}
  \item{bayes}{a logical indicating whether moderated estimates are used to construct the CD function.} 
  \item{p.nb}{a logical requesting a p-value based on a negative binomial model. If \code{TRUE}, the insertion-level p-values are derived from the DESeq (using the default setting),
  and then combined using the stouffer method.}
  \item{norm}{a logical indicating whether normalization is performed. If \code{TRUE} (default), TMM (trimmed mean of M values) normalization is used.} 
  \item{cut}{ an insertion is excluded from the analysis if the total read counts over all samples is less than cut. \code{cut}=0 (Default) include all data.} 
}

\value{
There are two output files, one is named as resTable, and the other is named as est.insertion.
 The resTable file is a data frame, which contains the following columns:
 \item{ID}{gene names}
 \item{NumInser}{number of insertions for each gene in input condition.}
 \item{Unique.NumInser}{number of unique insertions for each gene in input condition. If there is only one pool, Unique.NumInser should be equal to NumInser.}
 \item{Mean.x}{averaged counts for samples in condition x. The average is calculated after counts are normalized in DESeq.}
 \item{FC}{the fold change of condition B over A derived from the combined CD function (if replicates=FALSE, the FC is the ratio of the total counts within the gene).}
 \item{logFC}{log2 fold change of group B over group A.}
 \item{pvalue}{p-values indicating significances of the differential tests.}
  The est.insertion file is a data frame when there are replicates in all pools (for advanced users), and it contains estimates from the linear model for each insertion:
 \item{ID}{gene names}
 \item{Mean1,Mean2}{means of the two conditions.}
 \item{df}{degrees of freedom for the slope estimate.}
 \item{bhat}{slope estimate (logFC).}
 \item{se}{standard error of the slope estimate.}
}

\references{
Wang, H. and Song, M. (2001). Ckmeans.1d.dp: optimal k-means clustering in one dimension by
dyanamic programing. \emph{The R Journal}, 3(2), 29-33.


Zhao, L., Wu, W., Anderson, M. T., Li, Y. Mobley, H. L. T., and Bachman, M. A.(2017). 
Inseq: Identification of Conditionally Essential Genes in Transposon Sequencing Studies. 
Submitted to BMC Bioinformatics.
}


\author{
Lili Zhao \email{zhaolili@umich.edu}
}


%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{

%\dontrun{
data(serratia)

# Test the first 100 insertions
serr=serratia[1:100,]
# obtain the count matrix
countData=serr[,-c(1,2,3)]
condition=c(rep("Input",4),rep("Output",8))
pool=c(1,1,2,2,1,1,1,1,2,2,2,2)
geneID=as.character(serr$GeneID)
location=serr$Loc

foo<-TnseqDiff(countData, geneID, location, pool, condition)  

res=foo$resTable

# adjust pvalues using Benjamini & Hochberg
res$padj=p.adjust(res$pvalue, method = "BH")


## when no replicate in both conditions in each pool
countData=serr[,c(4,6,8,12)]
condition=c("Input","Input","Output","Output")
pool=c(1,2,1,2)
geneID=as.character(serr$GeneID)
location=serr$Loc

foo<-TnseqDiff(countData, geneID, location, pool, condition)  

res=foo$resTable


# when there is only one pool
countData=serr[,c(4,5,8,9,10,11)]
condition=c(rep("Input",2),rep("Output",4))
pool=c(1,1,1,1,1,1)
geneID=as.character(serr$GeneID)
location=serr$Loc

foo<-TnseqDiff(countData, geneID, location, pool, condition)  

res=foo$resTable

%}

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
