#'Summary estimate compiler
#'
#'Combining \code{\link{AlleleComp}} outputs from different samples with the variant sequence data.
#'
#'@param x A list object needs to be created by \code{\link{split}} from the sequencing data.
#'@param y A character vector of sample names or IDs.
#'
#'@details The variant sequence data needs to be split by sample names or IDs for \code{x}. And the input of \code{y}
#'has to be in the same order as that of the split object. See \code{example} for more details.
#'
#'@return A \code{dataframe} object identical to the original variant data with an additional column named \code{segment}
#'signifying the allelic make up of each variant in the corresponding sample.
#'
#'@seealso \code{\link{AlleleComp}}
#'
#'@examples
#' \donttest{#NB<-split(Neuroblastoma,Neuroblastoma$Sample)
#' #NB<-match.maker(x=NB,y=c("metastasis.1","metastasis.2","primary.1","primary.2"))
#' #View(NB)}
#' @export
match.maker<-function(x,y){
  dat<-data.frame()
  pb <- txtProgressBar(min = 0, max = 4, style = 3)
  for(n in 1 : length(x)){
    case<-x[[n]]
    match<-get(y[n])[[1]]
    match2<-get(y[n])[[2]]
    seg<-seg.track(match,match2)
    match$'allelic imbalance' <- seg[[1]]
    match$logR <- seg[[2]]
    case$logR <- case$'allelic imbalance' <- case$segment <- NA
    for(i in 1 : nrow(case)){
      for (j in 1 : nrow(match)){
        case$segment[i]<-ifelse(case$Chr[i] == match$chromosome[j] &&
                                  case$Pos[i] > match$start.pos[j] &&
                                  case$Pos[i] < match$end.pos[j], match$CN.profile[j], case$segment[i])

        case$'allelic imbalance'[i]<-ifelse(case$Chr[i] == match$chromosome[j] &&
                                              case$Pos[i] > match$start.pos[j] &&
                                              case$Pos[i] < match$end.pos[j], match$'allelic imbalance'[j], case$'allelic imbalance'[i])

        case$logR[i]<-ifelse(case$Chr[i] == match$chromosome[j] &&
                               case$Pos[i] > match$start.pos[j] &&
                               case$Pos[i] < match$end.pos[j], match$logR[j], case$logR[i])
      }
    }
    Sys.sleep(0.1)
    setTxtProgressBar(pb, n)

    dat<-rbind(dat,case)
    rm(list=c("case","match","match2"))
  }
  close(pb)
  return(dat)
}
