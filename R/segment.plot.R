#' Plot of allelic composition
#'
#' Departure in clusters of different allelic composition are portrayed for tumor sample 
#'
#' @param data A \code{match.maker} or \code{\link{AlleleComp}} derived object.
#' @param base.copy is the baseline  balanced copynumber present in the sample usually "1 + 1" or "2 + 2". 
#' 
#' @return A plot of the allelic segmentation with average log-transformed coverage ratios in X-axis and
#' average allelic-imbalances in the Y-axis. This plot can be interpreted in the similar fashion as described
#' by \href{https://genomebiology.biomedcentral.com/articles/10.1186/gb-2011-12-10-r108}{Rasmussen \emph{et al}., 2011}
#'
#' @seealso \code{\link{match.maker}}, \code{\link{AlleleComp}}
#' 
#' @examples \donttest{#segment.plot(data = data, base.copy = "1 + 1")}
#' 
#' @export

segment.plot<-function(data, base.copy){
  tmp1 = data
  tmp1$logR <- tmp1$logR - mean(subset(tmp1,tmp1$segment==base.copy)$logR)

  al<-vector()
  tmp1<-tmp1[order(tmp1$segment),]
  l<-unique(tmp1$segment)
  for (j in 1: length(l)){
    a1<-subset(tmp1,tmp1$segment==l[j])$logR
    if(length(a1) <2){
      al<-c(al,a1)
    } else {
      for (i in 1:length(a1)){
        mn<-median(a1[-i])
        al<-c(al,(mn-((mn-a1[i])/(length(a1)-1))))
      }
    }
  }
  tmp1$avg.logR<-al
  
  al<-vector()
  tmp1<-tmp1[order(tmp1$segment),]
  l<-unique(tmp1$segment)
  for (j in 1: length(l)){
    a1<-subset(tmp1,tmp1$segment==l[j])$'allelic_imbalance'
    if(length(a1) <2){
      al<-c(al,a1)
    } else {
      for (i in 1:length(a1)){
        mn<-median(a1[-i])
        al<-c(al,(mn-((mn-a1[i])/(length(a1)-1))))
      }
    }
  }
  tmp1$avg.Imb<-al

  g<-ggplot2::ggplot(tmp1,aes(x=avg.logR,y=avg.Imb))+
    geom_point(colour="grey", size=4, aes(alpha=0.5))+
    geom_point(aes(colour=factor(segment)), size=2)+
    labs(title="Scatter plot of VAF values",
         y="Allelic Imbalance",
         color="Allelic segments")+
    ylim(0,1)+
    xlim((min(tmp1$avg.logR) - 1.5),(max(tmp1$avg.logR) + 1.5))+
    theme(panel.background = element_rect(fill = "white"),
          plot.margin = margin(2, 2, 2, 2, "mm"),
          plot.background = element_rect(
            fill = "grey90",
            colour = "black",
            size = 1
          ))
  print(g)
}