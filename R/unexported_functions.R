Sys.setenv(`_R_S3_METHOD_REGISTRATION_NOTE_OVERWRITES_` = "0")
##================== Dot plot ==================
#Initial input and cluster.doc dot plot
dot_plot<-function(d){
  g<-ggplot(d,aes(x=sample,y=vaf))
  if("pred" %in% colnames(d)) {
    g<-g + geom_point(aes(color=factor(pred)),alpha=0.3, size=7)+
      geom_point(color='gray45', size=1)
  } else {
    g<-g + geom_point(aes(color=factor(sample)),alpha=0.3, size=7)+
      geom_point(color='gray45', size=1)+
      theme(legend.position="none")
    }
  g<-g+labs(title=NULL,
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
    scale_y_continuous(limits = c(0,1), expand = c(0, 0))
  return(g)
}
##================== Data extraction from VCF ==================
#Extracting allele frequencies from tumor and normal samples
#The vcfR file must ave the AD FORMAT input
#Adding chromosome, base pair position and reference allele, also
#performing general quality control steps
data.prep<-function(x, AD){
  ad <- as.data.frame(extract.gt(x,  element = AD)) #Approx read depth
  ad.normal.tmp<-unlist(strsplit(as.character(ad[,1]),split="[,]"))
  ad.tumor.tmp<-unlist(strsplit(as.character(ad[,2]),split="[,]"))
  ref.normal.tmp<-as.numeric(ad.normal.tmp[c(TRUE,FALSE)])
  alt.normal.tmp<-as.numeric(ad.normal.tmp[c(FALSE,TRUE)])
  ref.tumor.tmp<-as.numeric(ad.tumor.tmp[c(TRUE,FALSE)])
  alt.tumor.tmp<-as.numeric(ad.tumor.tmp[c(FALSE,TRUE)])
  seq<-data.frame(cbind(ref.normal.tmp,alt.normal.tmp,ref.tumor.tmp,alt.tumor.tmp))
  rm(list=ls()[grep("tmp",ls())])
  colnames(seq) <- c('AN','BN','AT','BT')
  seq$chromosome<-noquote(getCHROM(x))
  seq$position<-as.numeric(getPOS(x))
  seq$base.ref<-noquote(getREF(x))
  seq <- subset(seq,seq$BT>10 & seq$BN>10 & seq$AT>10 & seq$AN>10)
  return(seq)
}

##================== BAF and logR ==================
#Calcularing B allele frequency and relaive coverage
#for the entire genome per probe
allele.munge<-function(x){
  seq.genome<-data.frame()
  for (i in 1:22){
    seq.chr<-subset(x,x$chromosome==as.character(i))
    ##ratio of median of unscaled tumor and normal total coverage
    Rdep<-median(seq.chr$AT+seq.chr$BT)/median(seq.chr$AN+seq.chr$BN)
    for (j in 1:nrow(seq.chr)){
      ##B-allele frequency
      seq.chr$Bf[j]<-seq.chr$AT[j]/(seq.chr$AT[j]+seq.chr$BT[j])
      ##Scaled R
      seq.chr$R[j]<-(seq.chr$AT[j]+seq.chr$BT[j])/(seq.chr$AN[j]+seq.chr$BN[j])/Rdep
    }
    ##Creating a merged data
    seq.genome<-rbind(seq.genome,seq.chr)
    rm(seq.chr)
  }
  return(seq.genome)
}

##================== Rescale logR ==================
#Rescaling relative coverage with loess regression
#correcting on GC content of the sequence reads.
#This function creates two plots for each chromosome.
#One for original covergae and the other for rescaled.
#The plot outputs are sunk in a file specified by filename.
allele.summary<-function(x, filename){
  cat("correcting for GC bias",date(),"\n")
  file = paste0(filename,"_GCcorrection.pdf")
  pdf(file)
  seq.genome.fix<-data.frame()
  ##progress bar as it is done chromosome wise
  pb <- txtProgressBar(min = 0, max = 22, style = 3)
  for (i in 1:22){
    seq.chr<-subset(x,x$chromosome==as.character(i))
    ##GC bias correction
    tmp<-GC.correction(seq.chr ,i=i)
    ##Create the main data with corrected R
    seq.genome.fix<-rbind(seq.genome.fix,tmp)
    Sys.sleep(0.1)
    setTxtProgressBar(pb, i)
    rm(seq.chr)
  }
  dev.off()
  ##Creating corrected logR
  seq.genome.fix$logR<-suppressWarnings(log2(seq.genome.fix$R.fix))
  return(seq.genome.fix)
}

##================== CN estimation falconx ==================
#Copynumber estimates are generated with falcon
#Plots of B allele frequencies, relative coverage and
#allele-specific copy number estimates are sunk in
#pdf named with the filename iput
CN.summary<-function(x, filename, uniform.break){
  file = paste0(filename,"_CNprofile.pdf")
  #pdf(file)
  segments<-data.frame()
  cn_all<-list()
  for (i in 1:22){
    seq.chr<-subset(x,x$chromosome==as.character(i))
    seq.chr[seq.chr==0] <- NA
    seq.chr <- na.omit(seq.chr)
    seq.chr <- subset(seq.chr,seq.chr$BT>10 & seq.chr$BN>10 & seq.chr$AT>10 & seq.chr$AN>10)
    if (missing(uniform.break)){
      break.len = NULL } else {
        break.len = get.breakpoints(seq.chr$position,window=uniform.break)
      }
    cn <- falcon::getASCN(seq.chr[,c('AN','BN','AT','BT','R.fix')], tauhat=break.len)
    if (missing(uniform.break)){
      #view.falcon(cn,pos=as.numeric(seq.chr$position),i=i)
      tmp.pos<-c(min(seq.chr$position),seq.chr$position[cn$tauhat],max(seq.chr$position))
      chromosome<-seq.chr$chromosome[1:(length(cn$tauhat)+1)]
      start.pos<-(tmp.pos+1)[1:(length(cn$tau)+1)]
      end.pos<-tmp.pos[-1]
      segment.dat<-data.frame(chromosome,start.pos,end.pos)
      segment.dat$segment.length<-segment.dat$end.pos - segment.dat$start.pos
      tmp.cn<-cn$ascn
      segment.dat$CNt<-tmp.cn[1,]+tmp.cn[2,]
      segment.dat$A<-tmp.cn[1,]
      segment.dat$B<-tmp.cn[2,]
      segment.dat$CN.profile<-paste(ceiling(tmp.cn[2,]),"+",ceiling(tmp.cn[1,]))
      segments<-rbind(segments,segment.dat)
      rm(list=ls()[grep("tmp",ls())])
    } else {
      real.length<-length(unique(break.len))
      tmp.pos<-seq.chr$position[c(1,unique(break.len))]
      chromosome<-seq.chr$chromosome[1:real.length]
      start.pos<-(tmp.pos)[1:real.length]
      end.pos<-tmp.pos[-1]
      segment.dat<-data.frame(chromosome,start.pos,end.pos)
      segment.dat$segment.length<-segment.dat$end.pos - segment.dat$start.pos
      segment.dat$cns1<-cn$ascn[1,]
      segment.dat$cns2<-cn$ascn[2,]
      segments<-rbind(segments,segment.dat)
    }
    cn_all[[i]]<-cn
  }
  #dev.off()
  return(list(segments = segments, copy.num = cn_all))
}

#Centering the coverage ratio (R) on the unbalanced
#GC content
balancedCenter <- function(pos, numReads){
  center = 0
  netChange = sum((pos-center)*numReads)
  while(netChange > 0){
    center = center + 1
    netChange = sum((pos-center)*numReads)
  }

  #now zero in on zero (within 0.1%)
  margin = sum(numReads)*0.001
  adj = 0.5
  center = (center-1)+adj
  netChange = sum((pos-center)*numReads)
  while((abs(netChange) > margin) & (adj > 0)){
    adj = adj/2
    if(netChange > 0){
      center = center + adj
    }else{
      center = center - adj
    }
    netChange = sum((pos-center)*numReads)
  }
  return(center)
}

#Plotting the copy numbers
view.falcon <- function (output, pos = NULL, rdep = NULL, plot = "all", i, ...) {
  readMatrix = output$readMatrix
  tauhat = output$tauhat
  ascn = output$ascn
  AN = readMatrix$AN
  BN = readMatrix$BN
  AT = readMatrix$AT
  BT = readMatrix$BT
  Rf = as.numeric(readMatrix$R.fix)
  N = length(AT)
  tau = sort(unique(c(1, tauhat, N)))
  ascn1 = ascn[1, ]
  ascn2 = ascn[2, ]
  if (is.null(pos))
    pos = 1:N
  main.title = paste("Chromosome",i)
  myxlab = "Position (bp)"
  if (is.null(rdep))
    rdep = median(AT + BT)/median(AN + BN)
  if (plot == "all") {
    par(mfrow = c(3, 1))
  }
  if (plot == "all" || plot == "Afreq") {
    plot(pos, AN/(AN + BN), main = main.title, ylim = c(0, 1),
         ylab = "A freq", col = "gray", pch = ".",
         ...)
    points(pos, AT/(AT + BT), pch = ".", ...)
    abline(h = 0.5, col = "green")
    abline(v = pos[tau], col = "purple", lty = 4)
  }
  if (plot == "all" || plot == "RelativeCoverage") {
    #    if(is.null(Rf) == FALSE){
    plot(pos, Rf, ylab = "Relative Coverage",
         xlab = myxlab, pch = ".", ...)
    abline(h = 1, col = "green")
    abline(v = pos[tau], col = "purple", lty = 4)
    #    } else {
    #      plot(pos, (AT + BT)/(AN + BN)/rdep, ylab = "Relative Coverage",
    #           xlab = myxlab, pch = ".", ...)
    #      abline(h = 1, col = "green")
    #      abline(v = pos[tau], col = "purple", lty = 2)
    #    }
  }
  if (plot == "all" || plot == "ASCN") {
    falcon::view(output, pos = pos, plot="ASCN")
  }
}

#Getting 1Mb wide break points given a list of basepairs
#when uniform.break is meantioned in the wrapper function
get.breakpoints <- function(pos,window){
  l<-vector()
  tau<-c(seq(from=range(pos)[1], to=range(pos)[2], by=window), range(pos)[2])
  for (j in 2: length(tau)) {
    l<-c(l,which.min(abs(as.numeric(pos)-tau[j])))
  }
  return(l)
}

##================== GC content scaling ==================
GC.correction<- function(sequence, i){
  tumor<-sequence
  GC$pos<-paste0(GC$chr,":",GC$interval1,"-",GC$interval2)
  GC<-GC[,c(5,4)]
  window<-function(x){options(scipen=999);ceiling(x/10^4)*10^4}
  #Basepair positions has to be saved in the collumn name position
  interval.1<-window(tumor$position)
  #chromosomes has to be 1 to 22 integer saved in collumn name chromosome
  tumor$pos<-paste0('chr',tumor$chromosome,':',1+interval.1-10^4,"-",interval.1)
  tumor<-dplyr::left_join(tumor,GC,by='pos')
  tumor$GC.percent<-ceiling(tumor$GC*100)
  tumor<-subset(tumor, select=-c(pos,GC))

  #Performing loes correction on the relative coverage
  #which needs to be saved as the collumn name R
  loessMod30 <- loess(R ~ GC.percent, data=tumor, span=0.3)
  bmed=balancedCenter(sort(unique(loessMod30$fitted)),as.numeric(table(loessMod30$fitted)))
  reads.adj = loessMod30$fitted-bmed
  tumor$R.fix = tumor$R-reads.adj
  loessMod30.fix <- loess(R.fix ~ GC.percent, data=tumor, span=0.3)
  plot(y=tumor$R, x=tumor$GC.percent, type="p", main=paste("Loess Smoothing and Prediction for chromosome ", i),
       xlab="GC content percentage",
       ylab="relative coverage (R)",col="grey")
  points(loessMod30$fitted, x=tumor$GC.percent, pch= 16, col="blue")
  plot(y=tumor$R.fix, x=tumor$GC.percent, type="p",
       xlab="GC content percentage",
       ylab="Scaled relative coverage (R)",col="grey")
  points(loessMod30.fix$fitted, x=tumor$GC.percent, pch= 16, col="red")
  return(tumor)
}

##================== Wrapper for plot ==================
wrapper.plot <- function(plot.data, cn, estimated_copy, filename){
  #===================RELATIVE COVERAGE PLOT=====================
  g1 <- ggplot(plot.data,aes(x=position,y=R.fix))+
    geom_point(colour="royalblue4", size=1)+
    labs(title=NULL,
         subtitle = NULL,
         y="Relative coverage",
         x="Position (Mb)")+
    theme_light()+
    theme(legend.position="none")+
    ggpubr::grids(linetype = "dashed")+
    scale_x_discrete(limits=seq(0,1e9,by = 2e7), labels = seq(0,250, by = 20), expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0))
  g1 <- g1 + geom_hline(yintercept=1.0, color = "darkseagreen4", size=1)
  g1 <- g1 + theme(axis.title.x = element_text(size=14, face="bold"),
                   axis.title.y = element_text(size=14, face="bold"))

  #====================ALLELE FREQUENCY PLOT=======================
  plot.data$Norm.AF<-plot.data$AN/(plot.data$AN + plot.data$BN)
  plot.data$Tumor.AF<-plot.data$AT/(plot.data$AT + plot.data$BT)
  plot.data.new<-plot.data[,c('position','Norm.AF','Tumor.AF')]
  plot.data.new <- plot.data.new %>% gather("id", "value", 2:3)
  g2 <- ggplot(plot.data.new,aes(x=position,y=value))+
    geom_point(data = subset(plot.data.new,plot.data.new$id == 'Norm.AF'), colour= "grey", size=1)+
    geom_point(data = subset(plot.data.new,plot.data.new$id == 'Tumor.AF'), colour= "royalblue4", size=1)+
    labs(title=NULL,
         subtitle = NULL,
         y="Allele frequency",
         x="Position (Mb)")+
    theme_light()+
    theme(legend.position="none")+
    ggpubr::grids(linetype = "dashed")+
    scale_x_discrete(limits=seq(0,1e9,by = 2e7), labels = seq(0,250, by = 20), expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0,1))
  g2 <- g2 + geom_hline(yintercept=0.5, color = "darkseagreen4", size=1)
  g2 <- g2 + theme(axis.title.x = element_text(size=14, face="bold"),
                   axis.title.y = element_text(size=14, face="bold"))

  #====================COPY NUMBERPLOT===============================
  g3<-suppressWarnings(copy.num(cn, pos = plot.data$position))

  #====================ARRANGE THE COPY NUMBER PLOTS======================
  figure1 <- ggarrange(g1, g2, g3,
                       labels = NULL,
                       ncol = 1, nrow = 3)+
    theme(plot.margin = unit(c(0.3,0,0,0.2), "cm"))

  #====================ALLELIC IMBALANCE PLOT=======================
  chr = plot.data$chromosome[1]
  estimated_copy$annotation <- ifelse(estimated_copy$Chromosome == chr , 1, 2)
  m1<-round(min(estimated_copy$logr)-0.4,1)
  m2<-round(max(estimated_copy$logr)+0.4,1)
  g4 <- ggplot(estimated_copy,aes(x=logr,y=ai))+
    geom_point(colour="darkseagreen", size=2, aes(alpha=0.8))+
    geom_point(data = subset(estimated_copy,estimated_copy$Chromosome == chr), colour= "darkorchid3",
               aes(colour=factor(annotation)), size=1.5)+
    labs(title=NULL,
         subtitle = NULL,
         y="Allelic Imbalance",
         x="Average log-ratio")+
    theme_light()+
    theme(legend.position="none")+
    ggpubr::grids(linetype = "dashed")+
    scale_x_continuous(breaks = seq(m1,m2, by = 0.1), limits = c(m1,m2), expand = c(0, 0)) +
    scale_y_continuous(breaks = round(seq(0,1, by = 0.1),1), limits = c(0,1), expand = c(0, 0))+
    theme(plot.margin = unit(c(0.5,0.2,0.2,0.1), "cm"))
  g4 <- g4 + theme(axis.title.x = element_text(size=14, face="bold"),
                   axis.title.y = element_text(size=14, face="bold"))
  #====================ARRANGE ALL PLOTS TOGETHER=========================
  final.figure <- ggarrange(g4, figure1,
                            labels = NULL,
                            ncol = 2, nrow = 1)
  final.figure <- final.figure + labs(title=filename,subtitle = paste0("Chromosome ",chr))+
    theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 14, face = "bold"))

  #====================PRINT================
  print(final.figure)
}

##================== Copy number plot ==================
copy.num<-function (cn, pos = NULL){
  AT = cn$readMatrix$AT
  tauhat = cn$tauhat
  ascn = cn$ascn

  N = length(AT)
  tauhat = sort(unique(c(1, tauhat, N)))
  ascn1 = ascn[1, ]
  ascn2 = ascn[2, ]

  #Define labels and break points in terms of base pair
  if (is.null(pos)) {
    pos = 1:N
    if (is.null(xlab)) xlab = "SNP #"
  } else {
    tauhat = pos[tauhat]
    if (is.null(xlab)) xlab = "Position (bp)"
  }
  if (is.null(ylab)) ylab = "Allele-specific CN"

  pos1 = pos2 = rep(1, N)
  K = length(tauhat) - 1
  m = match(tauhat[1:K], pos)
  if (K > 1) {
    for (i in 1:(K - 1)) {
      pos1[m[i]:(m[i + 1] - 1)] = ascn1[i]
      pos2[m[i]:(m[i + 1] - 1)] = ascn2[i]
    }
  }
  pos1[m[K]:N] = ascn1[K]
  pos2[m[K]:N] = ascn2[K]
  p1 = l1 = p2 = l2 = c()
  if (K > 1) {
    for (i in 1:(K - 1)) {
      if (ascn1[i] > 1) {
        p1 = c(p1, m[i]:(m[i + 1] - 1))
      }
      else if (ascn1[i] < 1) {
        l1 = c(l1, m[i]:(m[i + 1] - 1))
      }
      if (ascn2[i] > 1) {
        p2 = c(p2, m[i]:(m[i + 1] - 1))
      }
      else if (ascn2[i] < 1) {
        l2 = c(l2, m[i]:(m[i + 1] - 1))
      }
    }
  }

  temp<-data.frame(cbind(pos,pos1,pos2))
  g <- ggplot(temp,aes(x=pos,y=value))+
    geom_point(colour = "darkseagreen4", aes (y = pos1), size = 1)+
    theme_light()+
    ggpubr::grids(linetype = "dashed")+
    scale_x_discrete(limits=seq(0,1e9,by = 2e7), labels = seq(0,250, by = 20), expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, max(c(ascn1,ascn2)) + 0.1)) +
    labs(title=NULL,
         subtitle = NULL,
         y="Copy number",
         x="Position (Mb)")
  g <- g + geom_point(colour = "darkseagreen4", aes(y = pos2), size = 1)

  l.1.0<-list()
  l.1.1<-list()
  l.1.2<-list()
  l.2.0<-list()
  l.2.1<-list()
  l.2.2<-list()
  l.3.0<-list()
  l.3.1<-list()
  l.3.2<-list()
  l.4.0<-list()
  l.4.1<-list()
  l.4.2<-list()
  for (i in 1:K) {
    if (i == K) {
      ids = m[K]:N
    } else {
      ids = m[i]:(m[i + 1] - 1)
    }
    nids = length(ids)

    #STEP1_1
    if (ascn1[i] > 1) {
      l.1.0[[i]]<-data.frame(cbind(pos = pos[ids], a1 = pos1[ids]))
      g <- g + geom_point(data = l.1.0[[i]], aes(x = pos, y = a1), size = 1, color = "indianred3")

      if (i > 1 && ascn1[i - 1] > 1) {
        a = ascn1[i - 1]
      } else {
        a = 1
      }
      l.1.1[[i]]<-data.frame(p1 = pos[ids[1]], p2 = pos[ids[1]], p3 = a, p4 = ascn1[i])
      g <- g + geom_segment(data = l.1.1[[i]], aes(x = p1, y = p3, xend = p2, yend = p4), size = 1, color = "indianred3")

      if (i < K - 1 && ascn1[i + 1] > 1) {
        a = ascn1[i + 1]
      } else {
        a = 1
      }
      l.1.2[[i]]<-data.frame(p1 = pos[ids[nids]], p2 = pos[ids[nids]], p3 = a, p4 = ascn1[i])
      g <- g + geom_segment(data = l.1.2[[i]], aes(x = p1, y = p3, xend = p2, yend = p4), size = 1, color = "indianred3")
    }

    #STEP1_2
    else if (ascn1[i] < 1) {
      l.2.0[[i]]<-data.frame(cbind(pos = pos[ids], a1 = pos1[ids]))
      g <- g + geom_point(data = l.2.0[[i]], aes(x = pos, y = a1), size = 1, color = "royalblue4")

      if (i > 1 && ascn1[i - 1] < 1) {
        a = ascn1[i - 1]
      } else {
        a = 1
      }
      l.2.1[[i]]<-data.frame(p1 = pos[ids[1]], p2 = pos[ids[1]], p3 = a, p4 = ascn1[i])
      g <- g + geom_segment(data = l.2.1[[i]], aes(x = p1, y = p3, xend = p2, yend = p4), size = 1, color = "royalblue4")

      if (i < K - 1 && ascn1[i + 1] < 1) {
        a = ascn1[i + 1]
      }
      else {
        a = 1
      }
      l.2.2[[i]]<-data.frame(p1 = pos[ids[nids]], p2 = pos[ids[nids]], p3 = a, p4 = ascn1[i])
      g <- g + geom_segment(data = l.2.2[[i]], aes(x = p1, y = p3, xend = p2, yend = p4), size = 1, color = "royalblue4")
    }

    #STEP2_1
    if (ascn2[i] > 1) {
      l.3.0[[i]]<-data.frame(cbind(pos = pos[ids], a1 = pos2[ids]))
      g <- g + geom_point(data = l.3.0[[i]], aes(x = pos, y = a1), size = 1, color = "indianred3")

      if (i > 1 && ascn2[i - 1] > 1) {
        a = ascn2[i - 1]
      } else {
        a = 1
      }
      l.3.1[[i]]<-data.frame(p1 = pos[ids[1]], p2 = pos[ids[1]], p3 = a, p4 = ascn2[i])
      g <- g + geom_segment(data = l.3.1[[i]], aes(x = p1, y = p3, xend = p2, yend = p4), size = 1, color = "indianred3")
      if (i < K - 1 && ascn2[i + 1] > 1) {
        a = ascn2[i + 1]
      }
      else {
        a = 1
      }
      l.3.2[[i]]<-data.frame(p1 = pos[ids[nids]], p2 = pos[ids[nids]], p3 = a, p4 = ascn2[i])
      g <- g + geom_segment(data = l.3.2[[i]], aes(x = p1, y = p3, xend = p2, yend = p4), size = 1, color = "indianred3")
    }

    #STEP2_2
    else if (ascn2[i] < 1) {
      l.4.0[[i]]<-data.frame(cbind(pos = pos[ids], a1 = pos2[ids]))
      g <- g + geom_point(data = l.4.0[[i]], aes(x = pos, y = a1), size = 1, color = "royalblue4")
      if (i > 1 && ascn2[i - 1] < 1) {
        a = ascn2[i - 1]
      } else {
        a = 1
      }
      l.4.1[[i]]<-data.frame(p1 = pos[ids[1]], p2 = pos[ids[1]], p3 = a, p4 = ascn2[i])
      g <- g + geom_segment(data = l.4.1[[i]], aes(x = p1, y = p3, xend = p2, yend = p4), size = 1, color = "royalblue4")
      if (i < K - 1 && ascn2[i + 1] < 1) {
        a = ascn2[i + 1]
      }
      else {
        a = 1
      }
      l.4.2[[i]]<-data.frame(p1 = pos[ids[nids]], p2 = pos[ids[nids]], p3 = a, p4 = ascn2[i])
      g <- g + geom_segment(data = l.4.2[[i]], aes(x = p1, y = p3, xend = p2, yend = p4), size = 1, color = "royalblue4")
    }
  }
  g <- g + geom_hline(yintercept=1.0, color = "darkseagreen4", size=1)
  g <- g + theme(axis.title.x = element_text(size=14, face="bold"),
                 axis.title.y = element_text(size=14, face="bold"))
}

##================== Regularization of Log ratios PREDICTED from NGS ==================
logR_scale_back <- function(data){
  data$ID<-1:nrow(data)
  data<-data[order(data$Value),]
  data$sign<-data$Value/abs(data$Value)
  colnames(data)[4]<-'old_value'
  data$Value<-scales::rescale(abs(data$old_value),c(0,2))
  data$Value<-data$Value * data$sign
  data<-data[order(data$ID),]
  data<-data[,-c(4:6)]
  return(data)
}

##================== Preping the 'AlleleComp' outputs for TAPS ==================
taps_prep <- function(segments,stats){
  alf<-stats[,c("chromosome","position","position","Bf")]
  colnames(alf)<-c("Chromosome","Start","End","Value")
  alf$Chromosome<-paste0("chr",alf$Chromosome)
  alf<-na.omit(alf)
  allele_freq<-alf

  Log2<-stats[,c("chromosome","position","position","logR")]
  colnames(Log2)<-c("Chromosome","Start","End","Value")
  Log2$Chromosome<-paste0("chr",Log2$Chromosome)
  Log2<-na.omit(Log2)
  Log2<-logR_scale_back(Log2)

  segments$copy.num<- ceiling(segments$A)+ceiling(segments$B)
  segments$chromosome<-as.numeric(segments$chromosome)
  # segments$chr<-paste0("chr",segments$chromosome)
  offset<-segments$chromosome - c(0,segments$chromosome[-nrow(segments)])
  # a<-vector()
  # for (i in 1: nrow(Log2)) {
  #   a <-c(a, subset(segments,
  #                   segments$chr == Log2$Chromosome[i] &
  #                     (segments$start.pos-offset) <= Log2$Start[i] &
  #                     segments$end.pos >= Log2$End[i])$copy.num)}
  # Log2$copy.num<-a
  segments$start.pos<-segments$start.pos-offset

  qc1<-segments[,c("chromosome","start.pos","end.pos","segment.length", "copy.num")]    #####
  colnames(qc1)<-c("chromosome","start","end","length", "copy.num")                     #####
  qc1$length<-round(qc1$length/1e6,3)
  for (i in 1: nrow(qc1)){
    suppressWarnings(qc1$snps[i]<-nrow(filter(stats, chromosome == qc1$chromosome[i]&
                                                position>qc1$start[i] &
                                                position<qc1$end[i])))
  }

  qc2 <- data.frame()
  for ( i in 1:nrow(qc1)){
    d <- data.frame(chromosome=numeric(),
                    start=numeric(),
                    end=numeric(),
                    copy.num=numeric(),                         #####
                    snps=integer())
    div.1 <- ceiling(qc1[i,]$snps/400)
    div.2 <- ceiling(qc1[i,]$length/4)
    div<-max(div.1,div.2)
    fxd.win <- floor(1e06*qc1[i,]$length/div)

    if (qc1[i,]$length < 4 | qc1[i,]$snps < 400) {
      d<-qc1[i,c("chromosome","start","end","copy.num","snps")] #####
    } else {
      for (j in 1 : div){
        d[j,]$chromosome <- qc1[i,]$chromosome
        d[j,]$start <- qc1[i,]$start + ((j-1)*fxd.win)
        d[j,]$end <- d$start[j] + fxd.win - 1
        d[j,]$copy.num <- qc1[i,]$copy.num                      #####
        d[j,]$snps <- nrow(filter(stats, chromosome == d[j,]$chromosome &
                                    position>=d$start[j] &
                                    position<=d$end[j]))
      }
    }
    qc2<-rbind(qc2,d)
    suppressWarnings(row.names(qc2)<-1:nrow(qc2))
  }
  colnames(qc2)[1]<-"chr"
  qc2$chr<-paste0("chr",qc2$chr)
  return(list(region = qc2,
              allele_freq = allele_freq,
              LogR = Log2))
}


##================== Creation of data structure for TAPS from summary stats ==================
stats_est <- function (regs, alf, Log2){
  ##progress bar as it is done chromosome wise
  pb <- txtProgressBar(min = 0, max = 22, style = 3)
  for (j in 1 : 22){
    chromosome<-paste0("chr",j)
    regs.sub<-subset(regs,regs$chr == chromosome)
    set<-data.frame()

    for (i in 0:(nrow(regs.sub)-1)) {
      start<-regs.sub[i+1,]$start
      end<-regs.sub[i+1,]$end

      #Create matchinf sets of alf Log2 data
      alf.sub<-subset(alf,alf$Chromosome == chromosome)
      alf.sub<-subset(alf.sub,alf.sub$Start >= start  & alf.sub$End <= end)
      Log2.sub<-subset(Log2,Log2$Chromosome == chromosome)
      Log2.sub<-subset(Log2.sub,Log2.sub$Start %in% alf.sub$Start)

      alf.sub<-alf.sub[!duplicated(alf.sub$Start),]
      Log2.sub<-Log2.sub[!duplicated(Log2.sub$Start),]

      alf.sub<-alf.sub[alf.sub$Start %in% Log2.sub$Start,]
      Log2.sub<-Log2.sub[Log2.sub$Start %in% alf.sub$Start,]

      if (nrow(Log2.sub) != nrow(alf.sub)) stop("Dimension mismatch while merging logR and allelic Imbalance")
      coords<-1:nrow(Log2.sub)

      #Quality control with PCA
      if (length(coords)>1) {
        a<-as.matrix(cbind(Log2.sub$Value, abs(alf.sub$Value-0.5)))
        U<-prcomp(a, scale. = TRUE, rank. = 10)$x
        pts<-apply(U, 2, function(x) which( abs(x - mean(x)) > (2 * sd(x)) ))

        if (length(pts) != 0){
          if (typeof(pts) == "integer") {
            if (length(pts) == 2){
              outliers<-unique(c(pts[1],pts[2]))
            } else {
              outliers<-unique(c(pts[,1],pts[,2]))
            }
          } else {
            outliers<-unique(c(as.numeric(pts$PC1),as.numeric(pts$PC2)))
          }
          coords.upd<-coords[-outliers]
        } else {
          coords.upd<-coords
        }

        #Deprecated PCA plots
        #theme_set(bigstatsr::theme_bigstatsr(0.8))
        #qplot(U[, 1], U[, 2]) + coord_equal()
        #theme_set(bigstatsr::theme_bigstatsr(0.8))
        #qplot(U[-outliers, 1], U[-outliers, 2]) + coord_equal()

        km<-kmeans(abs(alf.sub[coords.upd,]$Value-0.5),2)$center
        ai<-round(min(km[,1])/max(km[,1]),4)

        #Re-center logR based on AI and copy number
        if ("copy.num" %in% colnames(regs.sub)) {
          if (regs.sub$copy.num[i+1] > 2 & ai > 0.4){
            logr<-round(summary(Log2.sub[coords.upd,]$Value)[[6]],5)
          } else {
            logr<-round(median(Log2.sub[coords.upd,]$Value),5)
          }
        } else {
          logr<-round(median(Log2.sub[coords.upd,]$Value),5)
        }

        #quality control on AI based on no. of SNPs used
        ai<-ifelse(length(coords)<10, NA ,ai)
        ai<-ifelse(ai >=0.9, 1.0, ai)
        q<-cbind(j, start, end, logr, ai, length(coords), length(coords.upd))

      } else {
        #when coords = 1
        coords<-1:nrow(Log2.sub)
        logr<-round(median(Log2.sub[coords,]$Value),5)
        ai<-NA
        q<-cbind(j, start, end, logr, ai, length(coords), length(coords))
      }
      set<-rbind(set,q)
    }
    assign(chromosome,set)
    Sys.sleep(0.1)
    setTxtProgressBar(pb, j)
  }
  chrom<-rbind(chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,
               chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22)
  colnames(chrom)[c(1,6,7)]<-c("Chromosome","SNPs","Curated_SNPs")
  return(chrom)
}

##================== Plotting TAPS ==================
plot.taps <- function(data, sample){
  filename<-paste0(sample,".pdf")
  pdf(filename, width=25, height=8)
  for(chr in 1:22){
    data$annotation <- ifelse(data$Chromosome == chr , 1, 2)

    g <- ggplot(data,aes(x=logr,y=ai))+
      geom_point(colour="darkseagreen", size=2, aes(alpha=0.8))+
      geom_point(data = subset(data,data$Chromosome == chr), colour= "darkorchid3",
                 aes(colour=factor(annotation)), size=1.5)+
      labs(title=sample,
           subtitle = paste0("Chromosome ",chr),
           y="Allelic Imbalance",
           x="Average log-ratio")+
      theme_classic()+
      theme(legend.position="none")+
      ggpubr::grids(linetype = "dashed")+
      scale_x_continuous(breaks = round(seq(-2,2, by = 0.1),2), limits = c(-2,2), expand = c(0, 0)) +
      scale_y_continuous(breaks = round(seq(0,1, by = 0.1),1), limits = c(0,1), expand = c(0, 0))

    g <- g +  theme(plot.title = element_text(hjust = 0.5, color = "darkorange", size = 12, face = "bold"),
                    plot.subtitle = element_text(hjust = 0.5, color = "blue"))
    print(g)
  }
  dev.off()

}

##================== Wrapper function for stats_est & plot.taps ==================
taps_sorcery <- function (region , allele_freq, LogR, sample=NULL){
  est <- stats_est(regs = region, alf = allele_freq, Log2 = LogR)
  est <- na.exclude(est)
  #plot.taps(data = est, sample = sample)
  return(est)
}
##================== clonality determinant function ==================
clonality<-function(h1, h2, center){

tot.clone<-h1+h2
clone.diff<-1/tot.clone
all.center<-seq(1:(tot.clone-1))
expected.centers<-clone.diff*all.center

map<-vector()
for (i in 1:length(center)){
  dif<-abs(center[i]-expected.centers)
  dif<-cbind(dif,all.center)
  if(nrow(dif) > 1) {dif<-dif[order(dif[,1]),][1,]}
  map<-rbind(map,dif)
}
map<-cbind(map,center)
map<-map[order(map[,3], decreasing = TRUE),]

clonality<-'Clone'
for (i in 2:nrow(map)){
  if (map[i,2] == map[i-1,2]) pred <- 'Sub-clone'
  else pred <- 'Clone'
  clonality <- c(clonality,pred)
}
map<-cbind(data.frame(map),clonality)[,-1]
return(map)
}

##================== Initial input for cluster.doc ==================
Init.input <- function(x.1, al.comp = NULL){
  if (is.null(al.comp)) {
    h <- (readline("What is the suspected allelic composition of the variants: \n
    example : 1+1, 2+1 etc.\n"))
  } else {
    h <- al.comp
  }
  h1<-unlist(strsplit(h,"+")[1])[1]
  h2<-unlist(strsplit(h,"+")[1])[3]
  h1 <- as.numeric(h1)
  h2 <- as.numeric(h2)
  ##Check how the sample-wise distribution looks
  colnames(x.1)<-c("sample","vaf")
  g<-dot_plot(d = x.1)
  g<-g+theme(legend.position="none")+ scale_color_brewer(palette="Pastel1")
  cat("  Plotting the clonal smear of VAFs.
  Please guesstimate the following for all samples judging by the plot \n")
  print(g)
  r <- (readline("How many clonal VAF clouds do you think are present: \n"))
  r <- as.numeric(r)
  m <- (readline("How many sub-clonal VAF clouds do you think are present: \n"))
  m <- as.numeric(m)
  list(h1=h1, h2=h2, r=r, m=m)
}

##================== Plotting BIC for cluster.doc ==================
plot.bic<-function(tmp){
  if (is.null(dim(tmp)) == TRUE){
    tmp<-as.data.frame(tmp$profile)
    colnames(tmp)<-'SIN'
    cluster<-tmp$cluster<-as.numeric(rownames(tmp))

    ggplot(data=tmp, aes(x=cluster, y=SIN, color=SIN))+
      geom_point(size=4, show.legend = F)+
      theme_classic()+
      theme(legend.position = "bottom",
            legend.direction = "horizontal")+
      geom_line(size=1.5, show.legend = F)+
      labs(title="Changes in Smin as number of clusters increases",
           y="Smin statistics")+
      scale_x_continuous(breaks=cluster)
  }
  else{
    cluster<-as.numeric(rep(rownames(tmp),2))
    # BIC<-c(tmp[,1],tmp[,2])
    # group<-c(rep('Expectation',nrow(tmp)),rep('Variance',nrow(tmp)))BIC<-tmp[,1]
    ## INCLUSION OF AIC
    BIC<-tmp[,1]
    BIC<-c(BIC,BIC+4*(log(320)-4))
    group<-c(rep('BIC',nrow(tmp)),rep('AIC',nrow(tmp)))
    tmp<-as.data.frame(BIC)
    tmp$'Cluster Numbers'<-cluster
    tmp$metric<-group

    ggplot(data=tmp, aes(x=cluster, y=BIC, group=metric, color=metric, shape=metric, linetype=metric))+
      geom_point(size=4)+
      guides(alpha = FALSE)+
      scale_color_brewer(palette="Paired")+
      theme_classic()+
      theme(legend.position = "bottom",
            legend.direction = "horizontal")+
      geom_line(size=1.5, show.legend = F)+
      # labs(title="Changes in BIC as number of clusters increases",
      #      y="Bayesian Information Criteria (BIC)",
      #      metric="Evaluation Metric")+
      labs(title="Changes in AIC and BIC as number of clusters increases",
        y="Information estimate",
        metric="Evaluation Metric")+
      scale_x_continuous(breaks=cluster)
  }
}
