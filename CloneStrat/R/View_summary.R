#' Copynumber plots
#'
#' Summaries of copynumber estimations and allelic segmentations described in four plots
#'
#' @param data An \code{\link{AlleleComp}} object of the sequencing calls.
#' @param filename a \code{character} denominating sample name.
#'
#' @return A pdf file with plots on twenty two autosomes. The left panel of each plot shows allelic imbalance against average segmental log ratios.
#' In the right panel from top to bottom three plots describes chromosomal relative coverage, allele frequency and estimated copynumber.
#'
#' @examples \donttest{#met1<-AlleleComp(data = metastasis_1, AD = "AD", method = "naive", file = "TEMP");
#' #View_summary(data = met1, filename = "NB_met1")}
#'
#' @export

View_summary <- function(data, filename=NULL){
  segments<-data[[1]][[1]]
  stats<-data[[2]]
  taps_dat<-taps_prep(segments = segments, stats = stats)
  region<-subset(taps_dat$region,taps_dat$region$snps>9)
  est<-taps_sorcery(region = region , allele_freq = taps_dat$allele_freq, LogR = taps_dat$LogR)

  if(is.null(filename) == TRUE) {
    pdf.filename = paste0('Summary.plot','.pdf')
  } else {
    pdf.filename = paste0(filename,'.pdf')
  }
  pdf(pdf.filename, width=25, height=10)
  for(i in 1:22) {
    wrapper.plot(plot.data = subset(stats, stats$chromosome == i), cn = data[[1]][[2]][[i]], estimated_copy = est, filename = filename)
  }
  dev.off()
}
