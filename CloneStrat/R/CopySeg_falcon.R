##-------------------------------------------------
#'Copynumber estimation
#'
#'NGS probes are extracted rom a vcfR object, scaled and bias corrected to optimize estimatio
#'of allelelic composition. This function can handle only a combination of one tumor sample with
#'a matched normal sample. Analysis is performed using the package \code{\link{falcon}}
#'
#'@param data A \code{vcfR} object with one normal and one tumor sample. The \emph{AD} element
#'of the \emph{FORMAT} field is a manadatory input
#'@param AD a \code{character} denoting \emph{ID} for depth of the reference allele.
#'@param file.name A \code{character} string. this name will be used to save the scaled and unscaled
#'relative coverage plot along with the final copy number estimate plot in the working directory
#'@param uniform.break A numeric value signifying fixed length of the genomic window. Each window
#'is considered as distinct chromosomal segment with edges being the break points for copy number
#'estimation.
#'
#'@return A list of two data frames that is further used to obtain the allelic segmentation plot
#'
#'@details This function uses \code{\link{falcon}} to estimate allele specific copy number of all
#'sequeneced probes. Subsequently sliding window algorithm is used to generate chromosomal segments
#'with precicted distinct copynumbers. The relative coverages are sclaed with GC content of the binned
#'windows
#'@seealso \href{https://doi.org/10.1093/nar/gks001}{Benjamini \emph{et al}., 2012}
#'with a loess regression \code{\link{loess}}.
#'@export

CopySeg_falcon <- function(data, AD, file.name, uniform.break){
  seq<-data.prep(x=data, AD=AD)
  seq.genome<-allele.munge(x=seq)
  seq.genome.fix<-allele.summary(x=seq.genome, filename = file.name)
  if (missing(uniform.break)){
    segments_stats<-CN.summary(x=seq.genome.fix, filename = file.name)
  } else {
    segment_stats<-CN.summary(x=seq.genome.fix, filename = file.name, uniform.break = uniform.break)
  }
  return(list(segments_stats = segments_stats, seq.genome = seq.genome.fix))
}
