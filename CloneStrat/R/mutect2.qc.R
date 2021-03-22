#' Quality Control on Mutect2 output
#'
#' A quality control (QC) and transformation on the WES output from the Mutect2 variant caller. This re-organizes the data in a way that is
#' friendlier for using in \emph{CloneStrat}
#' @param WES A \code{dataframe} of the Mutect2 output
#' @param sample.name a \code{vector} of sample names or IDs
#' @return A transformed \code{dataframe} usable in \emph{CloneStrat} that represents data on each variant of each sample in rows
#' @examples \donttest{#res<-mutect2.qc(WES,sample.name)}
#' @export
mutect2.qc<-function(WES,sample.name){
  if (missing(sample.name) | typeof(sample.name)!='character') stop ("missing sample names")
  foo.1<-function(x){grep(x,colnames(WES))}
  req.col<-as.vector(sapply(sample.name,foo.1))
  if (length(req.col) %% 3 != 0) stop ("mismatch in number of columns. Only REF, ALT & AF allowed")
  west.1<-as_tibble(WES[,as.vector(sapply(sample.name,foo.1))])
  west.1[west.1 < 0.01] <- NA
  loop.len<-length(req.col)/3
  WES[,req.col]<-west.1
  colnames(WES)[req.col]<-rep(c("Mut_REF","Mut_ALT","VAF"),loop.len)
  west.3<-data.frame()
  for (i in 1 : loop.len){
    west.2<-na.omit(WES[,-c(req.col[-c((3*i-2):(3*i))])]) ##@i
    west.2$Sample<-sample.name[i] ##@i
    west.3<-rbind(west.3,west.2)
    west.3<-arrange.vars(west.3, c("Sample"=1))
    west.3<-arrange.vars(west.3, c("VAF"=2))
  }
  return(west.3)
}

#res<-mutect2.qc(WES,sample.name)
