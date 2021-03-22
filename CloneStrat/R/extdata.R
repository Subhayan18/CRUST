#' Random number generated WES data for eight hypothetical samples
#'
#' Data generated with varying random normal probabilities.
#' Ideal llelic composition is assumed resulting in two separate distinct clouds of clones and sub-clones.
#'
#' @docType data
#'
#' @usage data(test.dat)
#'
#' @format An object of class \code{"dataframe"}
#'
#' @return \code{sample} is column of IDs corresponding to 8 distinct samples.
#' @return \code{vaf} denotes the variant allele frequencies of each variant (see \code{annotation}).
#' @return \code{CCF} are the cancer cell fractions of each sample.
#' @return \code{annotation} indicates corresponding variants for which observations are notes in each row.
#' Variants can be shared among several samples as well as be private mutation.
#'
#' @keywords datasets
#'
#' @examples
#' \donttest{data(test.dat)
#' table(test.dat$CCF)
#' table(test.dat$annotation)
#' hist(test.dat$vaf)}
"test.dat"

#' Human neuroblastoma data
#'
#' Exome sequencing data of human neuroblastoma tumor samples available in public library.
#'
#' @docType data
#'
#' @usage data(Neuroblastoma)
#'
#' @format An object of class \code{"dataframe"}
#'
#' @return \code{Sample} is column of IDs corresponding to 4 samples (2x primary and 2x metastasis).
#' @return \code{VAF} denotes the variant allele frequencies.
#' @return \code{RefseqID} annotates each of the variants.
#'
#' @keywords datasets
#'
#'@seealso \code{\link{primary_1}} \code{\link{primary_2}} \code{\link{metastasis_1}} \code{\link{metastasis_2}}
#' @seealso \href{https://doi.org/10.1038/s41588-018-0131-y}{Karlsson \emph{et al}., 2018}
#'
#' @examples
#' \donttest{data(Neuroblastoma)}
"Neuroblastoma"

#'Human neuroblastoma tumor sample
#'
#'DNA sample collected from a primary tumor site (different than that of \emph{primary_2}) was sequenced.
#'This is a pre-processing vcfR file used for variant calling.
#'
#' @docType data
#'
#' @seealso \code{\link{Neuroblastoma}} \code{\link{primary_2}} \code{\link{metastasis_1}} \code{\link{metastasis_2}}
#' @seealso \href{https://doi.org/10.1038/s41588-018-0131-y}{Karlsson \emph{et al}., 2018}
"primary_1"

#'Human neuroblastoma tumor sample
#'
#'DNA sample collected from a primary tumor site (different than that of \emph{primary_1}) was sequenced.
#'This is a pre-processing vcfR file used for variant calling.
#'
#' @docType data
#'
#' @seealso \code{\link{Neuroblastoma}} \code{\link{primary_1}} \code{\link{metastasis_1}} \code{\link{metastasis_2}}
#' @seealso \href{https://doi.org/10.1038/s41588-018-0131-y}{Karlsson \emph{et al}., 2018}
"primary_2"


#'Human neuroblastoma tumor sample
#'
#'DNA sample collected from a metastatic site (different than that of \emph{primary_1}) was sequenced.
#'This is a pre-processing vcfR file used for variant calling.
#'
#' @docType data
#'
#' @seealso \code{\link{Neuroblastoma}} \code{\link{primary_1}} \code{\link{primary_2}} \code{\link{metastasis_2}}
#' @seealso \href{https://doi.org/10.1038/s41588-018-0131-y}{Karlsson \emph{et al}., 2018}
"metastasis_1"

#'Human neuroblastoma tumor sample
#'
#'DNA sample collected from a metastatic site (different than that of \emph{primary_1}) was sequenced.
#'This is a pre-processing vcfR file used for variant calling.
#'
#' @docType data
#'
#' @seealso \code{\link{Neuroblastoma}} \code{\link{primary_1}} \code{\link{primary_2}} \code{\link{metastasis_1}}
#' @seealso \href{https://doi.org/10.1038/s41588-018-0131-y}{Karlsson \emph{et al}., 2018}
"metastasis_2"
