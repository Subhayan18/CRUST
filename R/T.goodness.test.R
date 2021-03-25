#' Test of fit of clonal deconvolution
#'
#' A chi square test to assess the \emph{goodness of fit} of the clonal : sub-clonal clouds.
#' This test can be used to obtain outliers that do not fit into the proposed clonal deconvolution space.
#' @param x A \code{dataframe} with the first three columns in the specific order: sample name or ID of a variant, variant allele frquencies (VAF)
#' and cancer cell fraction (CCF)
#' @return A list of two objects. \emph{x} is same as the input \code{dataframe} with addede columns named \emph{expected VAF_}, \emph{chi_sq_}
#' and \emph{P value_} corresponding to each cloud of clone : Sub-clone combination. \emph{rej} is a subset of \code{x} containing variants
#' that fail the test for at least one cloud.
#' @return \emph{expected VAF_} represents estimated variant allele frequencies for a given cloud.
#' @return \emph{chi_sq_} is the Chi square test statistic for the cloud.
#' @return \emph{P value_} is the P value corresponding to the \emph{chi_sq_} statistic.
#' @examples \donttest{
#' #test<-T.goodness.test(test.dat)
#' #head(test)}
#' @export
T.goodness.test<-function (x){
  h <- (readline("What is the suspected chromosomal segmentation profile of the sample: \n
  (WARNING! Include spaces between numbers and '+')\n
  example : 1 + 2, 2 + 2, 19 + 57 etc.\n"))
  h1<-unlist(strsplit(h," ")[1])[1]
  h2<-unlist(strsplit(h," ")[1])[3]
  h1 <- as.numeric(h1)
  h2 <- as.numeric(h2)
  colnames(x)[1:3]<-c("Sample","VAF","CCF")
  a = x[,paste0('expected VAF (',h2,' mutated allele)')] <- round((x$CCF*(h2))/(x$CCF*(h1+h2)+(1-x$CCF)*2), 2)
  a = x[,paste0('chi_sq (',h2,' mutated allele)')] <- round(2*x$VAF/a, 2)
  x[,paste0('P value (',h2,' mutated allele)')] <- round(pchisq(a,1,lower.tail = F), 3)
  for (i in 1 : h2-1) {
    a = x[,paste0('expected VAF (',h2-i,' mutated allele)')] <- round((x$CCF*(h2-i))/(x$CCF*(h1+h2)+(1-x$CCF)*2), 2)
    a = x[,paste0('chi_sq (',h2-i,' mutated allele)')] <- round(2*x$VAF/a, 2)
    x[,paste0('P value (',h2-i,' mutated allele)')] <- round(pchisq(a,1,lower.tail = F), 3)
  }
  suppressWarnings({rej = filter_all(as_tibble(x), any_vars(. < 0.05))})
  list(x=x,rej=rej)
}
