\name{serratia}
\alias{serratia}
\docType{data}
\title{
serratia data example
}
\description{
a small subset of the real serratia transposon sequencing data. serratia is a data frame containing insertion count data from two inoculum pools, and each pool has 6 samples
(2 input and 4 output samples).}
\usage{data("serratia")}
\format{
  A data frame with 2442 observations on the following 15 variables.
\describe{
    \item{GeneID}{gene name}
    \item{Length}{gene length}
    \item{Loc}{insertion locations}
    \item{SM_control_lib1_rep1}{first input sample from pool 1}
    \item{SM_control_lib1_rep2}{second input sample from pool 1}
    \item{SM_control_lib2_rep1}{first input sample from pool 2}
    \item{SM_control_lib2_rep2}{second input sample from pool 2}
    \item{SM_target_lib1_rep1}{first output sample from pool 1}
    \item{SM_target_lib1_rep2}{second output sample from pool 1}
    \item{SM_target_lib1_rep3}{third output sample from pool 1}
    \item{SM_target_lib1_rep4}{forth output sample from pool 1}
    \item{SM_target_lib2_rep1}{first output sample from pool 2}
    \item{SM_target_lib2_rep2}{second output sample from pool 2}
    \item{SM_target_lib2_rep3}{third output sample from pool 2}
    \item{SM_target_lib2_rep4}{forth output sample from pool 2}
  }
}

\examples{
data(serratia)

}

