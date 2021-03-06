% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/match.maker.R
\name{match.maker}
\alias{match.maker}
\title{Summary estimate compiler}
\usage{
match.maker(x, y)
}
\arguments{
\item{x}{A list object needs to be created by \code{\link{split}} from the sequencing data.}

\item{y}{A character vector of sample names or IDs.}
}
\value{
A \code{dataframe} object identical to the original variant data with an additional column named \code{segment}
signifying the allelic make up of each variant in the corresponding sample.
}
\description{
Combining \code{\link{AlleleComp}} outputs from different samples with the variant sequence data.
}
\details{
The variant sequence data needs to be split by sample names or IDs for \code{x}. And the input of \code{y}
has to be in the same order as that of the split object. See \code{example} for more details.
}
\examples{
\donttest{#NB<-split(Neuroblastoma,Neuroblastoma$Sample)
#NB<-match.maker(x=NB,y=c("metastasis.1","metastasis.2","primary.1","primary.2"))
#View(NB)}
}
\seealso{
\code{\link{AlleleComp}}
}
