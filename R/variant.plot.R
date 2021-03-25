#' Multi-sample variant plot
#'
#' Plotting a specific variant present in more than one WES sample
#' @param CD.obj A \code{cluster.doc object}
#' @param annotation.col name of the column containing annotations of the variants
#' in original WES \code{dataframe} used in the clonal deconvolution using \code{cluster.doc}
#' @param variant a \code{character} string specifying \emph{only one} annotation which is to be displayed
#' @return A plot object with the relevant annotation highlighted
#' @examples
#' \donttest{#cd.res<-cluster.doc(test.dat,1,2)
#' #variant.plot(cd.res,'annotation','variant_74')}
#' @export
variant.plot<-function(CD.obj,annotation.col,variant){
  if (missing(CD.obj) | typeof(CD.obj)!='list') stop ("missing CD.obj")
  if (missing(annotation.col) | typeof(annotation.col)!='character') stop ("missing annotation.col")
  if (missing(variant) | typeof(annotation.col)!='character') stop ("missing variant")
  x<-CD.obj$predicted.data
  x$marker<-ifelse(x[annotation.col] == variant , variant, as.character(x$'Predicted clonal structure'))
  x$marker<-factor(x$marker, levels=c(variant,'Clone','Sub-clone'))
  names<-CD.obj$'user input'
  ggplot(x,aes_string(x=names[1],y=names[2]))+
    geom_point(aes(color=factor(marker)),size = ifelse(x$marker == variant,7,6),
               alpha=ifelse(x$marker == variant,0.6,0.3))+
    geom_point(color='gray45', size=1)+
    labs(title="Scatter plot of VAF values",
         y="Variant allele frequency",
         color="Annotations")+
    scale_color_manual(values=c('aquamarine4',"tan1","mediumpurple3"))+
    labs(title=NULL,
    y="Variant allele frequency",
    color="samples")+
    theme_classic()+
    theme(axis.ticks.length = unit(0.15, "cm"))+
    theme(axis.title.y = element_text(size = 16),
          axis.text.x = element_text(face="bold", size=10, angle=30, hjust = 0.5, vjust = 0.5),
          axis.text.y = element_text(face="bold", size=12))+
    theme(legend.position = "bottom",
          legend.direction = "horizontal")+
    labs(y="Variant allele frequency",
         x=element_blank())+
    scale_y_continuous(limits = c(0,max(x[names[2]])+0.1), expand = c(0, 0))
}
