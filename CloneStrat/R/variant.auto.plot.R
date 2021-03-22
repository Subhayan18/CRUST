#' Automated Multi-sample plot
#'
#' Automated plotting of all variants present in the WES data
#' @param CD.obj A \code{cluster.doc object}
#' @param annotation.col name of the column containing annotations of the variants
#' in original WES \code{dataframe} used in the clonal deconvolution using \code{cluster.doc}
#' @return Plot objects with the relevant annotation highlighted.
#' @return This function plots all variants present in the sample. Depending on the number of variants this can generate a \emph{lot}
#' of plots. All of these plots will be saved under a new directory named \code{img} inside the working directory.
#' Hence, it is important to check that there are no directory named \code{img} inside the working directory
#' @examples
#' \donttest{#cd.res<-cluster.doc(test.dat,1,2)
#' #variant.auto.plot(cd.res,'annotation')}
#' @export
variant.auto.plot<-function(CD.obj,annotation.col){
  if (missing(CD.obj) | typeof(CD.obj)!='list') stop ("missing CD.obj")
  if (missing(annotation.col) | typeof(annotation.col)!='character') stop ("missing annotation.col")
  x<-CD.obj$predicted.data

  if (!dir.exists(file.path(getwd(),"img"))) {
    dir.create("./img") }
    else {
      stop("Directory already exists: ")
    }

  list.varied<-unique(x[annotation.col])
  for (i in 1 : nrow(list.varied)){
    variant=as.character(list.varied[i,1])
    x$marker<-ifelse(x[annotation.col] == variant , variant, x$'Predicted clonal structure')
    x$marker<-factor(x$marker, levels=c(variant,'Clone','Sub-clone'))
    names<-CD.obj$'user input'
    ggplot(x,aes_string(x=names[1],y=names[2]))+
      geom_point(aes(color=factor(marker)),size = ifelse(x$marker == variant,4,1.5),
                 shape=ifelse(x$marker == variant,17,16))+
      labs(title="Scatter plot of VAF values",
           y="Variant allele frequency",
           color="Annotations")+
      ylim(0,max(x[names[2]])+0.1)+
      scale_color_manual(values=c('palegreen2',"tan1","mediumpurple3"))+
      theme_light()
    ggsave(paste0(getwd(),"/img/",list.varied[i,1],".tiff"), units="in", width=8, height=6, dpi=600, compression = 'lzw')
  }
}
