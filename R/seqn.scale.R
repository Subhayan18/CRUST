#' Probabilistic quotient normalization of DNA sequencing data
#'
#' A normalization technique based on cancer / tumor cell fractions of the samples sequenced to infer homogeneity
#'
#' @param x A \code{dataframe} containing summary from DNA sequencing with first column as sample IDs of corresponding variants.
#' @param vaf The column number of \code{x} that includes VAFs.
#' @param CCF The column number of \code{x} that includes CCFs.
#'
#' @details Probabilistic quotient normalization normalization technique described in
#' \emph{Dieterle, et al. (2006)} applied on the cancer cell fraction (CCF)
#' of respective samples to rescale variant allele frequencies (VAF) accordingly. The general idea is to put most
#' confidence in the sample with highest CCF and adjust the VAFs of other samples based on the departure in CCF
#' of the other samples from that with the highest.
#' @details This method is particularly suggested if the CCFs accross samples vary more than 10%.
#'
#' @return A \code{dataframe} with all the elements of \code{x} with the new estimated VAFs in the column
#' \emph{scaled.vaf} and an additional column \emph{unscaled.vaf} that includes the original VAFs.
#'
#' @seealso \code{\link{cluster.doc}}
#' @examples
#' \donttest{#pqn.dat<-seqn.scale(test.dat,vaf=2,CCF=3)
#' #hist(pqn.dat$scaled.vaf)}
#' @export
seqn.scale <- function(x,vaf,CCF) {
  if(missing(vaf) | missing(CCF)) {
    stop("missing arguments!")
  }
  x$ID<-1:nrow(x)
  train<-x[x[[CCF]] == max(x[[CCF]]),]
  test<-dplyr::anti_join(x,train)
  N.pqn<-suppressWarnings(pqn(sample=test[,c(vaf,CCF)],ref=train[,c(vaf,CCF)]))
  X.pqn<-rbind(train,test)
  X.pqn$pqn_scaled_vaf<-c(train[[vaf]],N.pqn[,1])
  for (i in 1 : length(X.pqn$pqn_scaled_vaf)) {
    X.pqn$pqn_scaled_vaf[i] <- ifelse(X.pqn$pqn_scaled_vaf[i] >
                                        1, 2-X.pqn$pqn_scaled_vaf[i], X.pqn$pqn_scaled_vaf[i])
  }
  X.pqn$pqn_scaled_vaf<-round(X.pqn$pqn_scaled_vaf,4)
  X.pqn<-X.pqn[order(X.pqn$ID),]
  X.pqn<-subset(X.pqn, select=-c(ID))
  colnames(X.pqn)[vaf]<-'unscaled.vaf'
  X.pqn<-arrange.vars(X.pqn, c("pqn_scaled_vaf"=vaf+1))
  colnames(X.pqn)[vaf+1]<-'scaled.vaf'
  return(X.pqn)
}
