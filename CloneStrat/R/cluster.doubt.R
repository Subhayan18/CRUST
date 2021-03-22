#' User overriden clonal deconvolution
#'
#' Sample specific user curated Clone / Sub-clone decomposition of DNA sequencing data
#' @param CD.obj A \code{cluster.doc} object
#' @param sample Collumn number of the \emph{predicted.data} from the \code{cluster.doc}
#' output that contain sample IDs
#'@param vaf Collumn number of the \emph{predicted.data} from the \code{cluster.doc}
#' output that contain variant allele frequencies used for the analysis.
#' @param sample.name a vector of sample IDs
#' @param cluster.num a numeric vetor of clone/sub-clonal split of respective sample
#' @return A list of 3 objects
#' @return \code{fitted.cluster} includes the clustering results from the final fit with user input
#' @return \code{predicted.data} A dataframe shows the changed clustering results due to the user defined
#' clone / sub-clone smear for the selected samples
#' @seealso \code{\link{cluster.doc}}
#' @examples
#' \donttest{#cd.res<-cluster.doc(test.dat)
#' #cd.new<-cluster.doubt(cd.res,sample,vaf,c("Sample_1","Sample_3"),c(2,2,3,2))}
#' @export
cluster.doubt<-function(CD.obj,sample,vaf,sample.name,cluster.num){
  if (missing(CD.obj) | typeof(CD.obj)!='list') stop ("missing CD.obj")
  if (missing(sample.name)) stop ("missing sample.name")
  if (missing(cluster.num) | typeof(cluster.num)!='double') stop ("missing cluster.num")
  if (length(cluster.num) != 2*length(sample.name)) stop ("mismatch in number of clones and subclones")

  cluster.num.0<-cluster.num
  cluster.num.1<-unname(tapply(cluster.num, (seq_along(cluster.num)-1) %/% 2, sum))
  cluster.num.2<-unname(tapply(cluster.num, (seq_along(cluster.num)-1) %/% 2, paste0))
  for(i in 1 : length(cluster.num.2)) {cluster.num.2[i] <- paste0(cluster.num.2[[i]][1],"+",cluster.num.2[[i]][2])}
  cluster.num.2<-unlist(cluster.num.2)

  temp.save<-colnames(CD.obj$predicted.data)[c(sample,vaf)]
  colnames(CD.obj$predicted.data)[c(sample,vaf)]<-c('sample','vaf')
  dat<-CD.obj$predicted.data[,c(sample,vaf)]
  improv.0<-as.data.frame(cbind(sample.name,cluster.num.2))
  old.data<-suppressWarnings(anti_join(CD.obj$predicted.data,improv.0,by=c('sample'='sample.name')))
  old.data<-old.data[,c(sample,vaf,ncol(old.data))]
  improv.1<-split(improv.0, improv.0$cluster.num.2)
  improv.data<-data.frame()
  for (i in 1 : length(cluster.num.1)){
    improv.2<-improv.1[[as.character(cluster.num.2[i])]]
    data<-suppressWarnings(inner_join(dat,improv.2,by=c('sample'='sample.name'))[,c(1,2)])
    kmeansM<-kmeans(as.matrix(data[,2]),centers=cluster.num.1[i])
    tmp.1<-rank(kmeansM$centers)
    tmp.2<-as.numeric(tmp.1[as.factor(kmeansM$cluster)])
    max.tmp.2<-max(tmp.2)
    min.tmp.2<-min(tmp.2)
    if (cluster.num.1[i] %% 2 == 0) {
      data$'Predicted clonal structure'<- noquote(unlist(lapply(tmp.2 %% 2 == 0, ifelse, "Clone", "Sub-clone")))
    } else if (cluster.num.1[i] %% 2 != 0 & cluster.num[2*i] <= cluster.num[(2*i)-1]){
      data$'Predicted clonal structure'<- noquote(unlist(lapply(tmp.2 %% 2 == 0 | tmp.2 == max.tmp.2, ifelse, "Clone", "Sub-clone")))
    } else {
      data$'Predicted clonal structure'<- noquote(unlist(lapply(tmp.2 %% 2 != 0, ifelse, "Clone", "Sub-clone")))
      data$'Predicted clonal structure'<- ifelse(tmp.2 == min.tmp.2, "Sub-clone", data$'Predicted clonal structure')
    }
    improv.data=rbind(data,improv.data)
  }
  improv.data=improv.data[1:(nrow(improv.data)/length(cluster.num.1)),]
  data<-rbind(old.data,improv.data)

  #Print the clones and sub-clones
  colnames(data)<-c("sample","vaf","pred")
  data<-data[order(data$sample),]

  g <- dot_plot(data)
  g <- g + scale_color_manual(values=c("tan1","mediumpurple3"))
  print(g)


  CD.obj$predicted.data$'Predicted clonal structure' <- data$'pred'
  colnames(CD.obj$predicted.data)[c(sample,vaf)]<-temp.save

  invisible(list(fitted.cluster=kmeansM,
                 predicted.data = CD.obj$predicted.data,
                 'user input' = CD.obj$'user input'))
}
