#' Copynumber estimation
#'
#' Allelic segmentations are estimated for one sample at a time with unfiltered sequencing calls.
#'
#' @param data A \code{vcfR} object of the sequencing calls.
#' @param AD a \code{character} deoning \emph{ID} for depth of the reference allele.
#' This is often separately present in the VCF file. Default is \code{NULL}.
#' @param file.name an optional \code{character} to define output file name. Default is \emph{tumor.sample}.
#' @param method Algorithm to be used for copy number calculations. options include "apriori" wich
#' uses \code{\link{CopySeg_sequenza}} and "naive" using \code{\link{CopySeg_falcon}}.
#' @param uniform.break A numeric value signifying fixed length of the genomic window. Each window
#' is considered as distinct chromosomal segment with edges being the break points for copy number
#' estimation. A good window length is 1Mb (i.e. 1e6)
#'
#' @return A transformed \code{dataframe} usable in \emph{CloneStrat} that represents data on all variants
#' in the .vcf file. It returns summaries on the variants with the collumn \emph{CN.profile} depicting
#' the estimated allelic compositions.
#'
#' @details The function writes a \emph{.txt} data in working directory with the name defined in \code{file.name} used by \emph{sequenza}.
#' The output file written can be used in conjunction with post variants call sequence file. These can be merged and used for surther analysis
#' with \code{\link{cluster.doc}} or \code{\link{seqn.scale}}
#'
#' @seealso \code{\link{segment.plot}}
#'
#' @examples
#' \donttest{#AlleleComp(data = x, AD = "AD", method = "naive")}
#'
#' @export

AlleleComp<-function(data, AD, file.name, method, uniform.break){
  if (missing(data) | class(data)[1] != 'vcfR')
    stop ("missing or incorrect entry for x")
  if (missing(AD) | class(AD) != 'character')
    stop ("missing or incorrect entry for AD")
  if (!missing(file.name) & class(file.name) != 'character')
    stop ("incorrect entry for file.name")
  if (missing(file.name)) {file.name="test"}
  if (missing(method) & method != "apriori" & method != "naive" & class(method) != 'character')
    stop ("A valid entry for method is required")

  if (method == 'apriori'){
    out<-CopySeg_sequenza(x=data, AD=AD, file.name=file.name)
  } else if (method == 'naive'){
    if (missing(uniform.break)){
      out<-CopySeg_falcon(data=data, AD=AD, file.name=file.name)
    } else{
      out<-CopySeg_falcon(data=data, AD=AD, file.name=file.name, uniform.break=uniform.break)
    }
  } else stop("Invalid input provided for method")
  return(out)
}
