#' Copynumber estimation
#'
#' Allelic segmentations are estimated for one sample at a time with unfiltered sequencing calls
#' using the package \code{sequenza}.
#' This function can handle only a combination of one tumor sample with
#'a matched normal sample.
#'
#' @param x A \code{vcfR} object of the sequencing calls.
#' The sample names can be queried from \code{x}.
#' @param AD a \code{character} deoning \emph{ID} for depth of the reference allele.
#' This is often separately present in the VCF file. Default is \code{NULL}.
#' @param file.name an optional \code{character} to define output file name. Default is \emph{tumor.sample}.
#'
#' @return A transformed \code{dataframe} usable in \emph{CloneStrat} that represents data on all variants
#' in the .vcf file. It returns summaries
#' on the variants with the collumn \emph{CN.profile} depicting the estimated allelic segmentations.
#'
#' @details The function writes a \emph{.txt} data in working directory with the name defined in
#' \code{file.name} used by \emph{sequenza}. The output file written can be used in conjunction with
#' post variants call sequence file. These can be merged and used for surther analysis
#' with \code{\link{cluster.doc}} or \code{\link{seqn.scale}}

CopySeg_sequenza<-function(x, AD, file.name){
  ad <-as.data.frame(extract.gt(x,  element = AD, as.numeric = FALSE))
  ad.normal.tmp<-unlist(strsplit(as.character(ad[,1]),split="[,]"))
  ad.tumor.tmp<-unlist(strsplit(as.character(ad[,2]),split="[,]"))
  ref.normal.tmp<-as.numeric(ad.normal.tmp[c(TRUE,FALSE)])
  alt.normal.tmp<-as.numeric(ad.normal.tmp[c(FALSE,TRUE)])
  ref.tumor.tmp<-as.numeric(ad.tumor.tmp[c(TRUE,FALSE)])
  alt.tumor.tmp<-as.numeric(ad.tumor.tmp[c(FALSE,TRUE)])

  chromosome<-noquote(getCHROM(x))
  position<-as.numeric(getPOS(x))
  base.ref<-noquote(getREF(x))
  depth.normal=alt.normal.tmp+ref.normal.tmp
  depth.tumor=alt.tumor.tmp+ref.tumor.tmp
  depth.ratio<-round((alt.tumor.tmp+ref.tumor.tmp)/(alt.normal.tmp+ref.normal.tmp),3)
  Bf<-round(alt.tumor.tmp/(alt.tumor.tmp+ref.tumor.tmp),3)
  Af<-1-Bf
  zygosity.normal<-ifelse(Bf < 0.9 & Bf > 0.1, 'het', 'hom')

  tumor<-dplyr::as_tibble(cbind(chromosome, position, base.ref, depth.normal,
                                depth.tumor, depth.ratio, Af, Bf, zygosity.normal))
  GC$pos<-paste0(GC$chr,":",GC$interval1,"-",GC$interval2)
  GC<-GC[,c(5,4)]
  window<-function(x){options(scipen=999);ceiling(x/10^4)*10^4}
  interval.1<-window(position)
  tumor$pos<-paste0('chr',chromosome,':',1+interval.1-10^4,"-",interval.1)
  tumor<-dplyr::left_join(tumor,GC,by='pos')

  tumor$GC.percent<-ceiling(tumor$GC*100)
  tumor<-subset(tumor, select=-c(pos,GC))

  tumor$good.reads<-ceiling(0.95*as.numeric(tumor$depth.tumor))
  tumor$AB.normal<-'.'
  tumor$AB.tumor<-'.'
  tumor$tumor.strand<-'.'
  tumor<-tumor[tumor$chromosome %in% as.character(1:22),]
  tumor<-subset(tumor,tumor$depth.normal != '0')
  tumor<-na.omit(tumor)

  name.txt<-paste0("sqz_",file.name,".txt")

  print("Writing the sequenza data file")
  write.table(tumor,name.txt,row.names=F,quote=F,sep="\t")
  if (nrow(tumor) > 10000) {print("Adequate variants present")}
  else {stop("Inadequate (less than 10k) variants present in the file! Please make sure the constitutional reads are present.")}

  test<-sequenza.extract(name.txt, verbose=FALSE)

  print("Evaluating summary statistics")
  CP <- sequenza.fit(test)
  print("Evaluating mode of allelic composition")
  confint <- get.ci(CP)
  seqz.data <- read.seqz(name.txt)
  test.CN<-as_tibble(baf.bayes(seqz.data$Bf, seqz.data$depth.ratio,
                               confint$max.cellularity, confint$max.ploidy,avg.depth.ratio = 1))
  tumor$CN.profile<-paste(test.CN$A,"+",test.CN$B)
  tumor$CN.profile<-ifelse(tumor$CN.profile == "0 + 0",NA,tumor$CN.profile)
  tumor<-tumor[,-c(12,13,14)]

  test.CN <- ChromSeg(sequenza.extract = test,cp.table = CP)
  test.CN$CN.profile<-paste(test.CN$A,"+",test.CN$B)
  test.CN<- na.omit(subset(test.CN,test.CN$CN.profile != "0 + 0"))
  print("Done")

  list(test.CN,tumor)
}
