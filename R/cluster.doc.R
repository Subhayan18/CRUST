#' Clonal deconvolution
#'
#' Clone / Sub-clone decomposition of DNA sequencing data. This is recommended to be used for more than one sample preferably collected from
#' the same individual at different times. If the sample qualities vary, it is recommended to perform scaling first with \code{\link{seqn.scale}}.
#'
#' @param data A \code{dataframe} containing summary from DNA sequencing. It must include a column of sample IDs
#' and a corresponding column with the variant allele frequencies.
#' @param sample \code{Integer or character} of the column name or column number of the sample IDs.
#' @param vaf \code{Integer or character} of the column name or column number of the variant allele frequency.
#' @param allele.comp \code{Character} string for allelic composition of the variants. example: '1+1' or '2+3' etc.
#' @param n.clone Optional \code{integer} for number of suspected clones, default NULL.
#' @param n.subclone Optional \code{interger} for number of suspected subclones, default NULL.
#' @param optimization.method Method to find optimal number of clusters; \emph{GMM} or \emph{bootstrap}. Default is \emph{GMM}.
#' @param clustering.method Clustering methods; \emph{HKM}, \emph{bootkm} or \emph{hybrid}. Default is \emph{hkm}.
#' @param clonality Method for determining clonality of the predicted clusters; \emph{Allelic composition} (default) or \emph{density}
#' @param instruct \code{Character} input for accepting program suggestion.
#'
#' @details \code{cluster.doc} is meant to do two things, first determine the optimum number of clusters that \emph{should} be fitted and
#' second, to infer what groups the clusters thus obtained should be assigned to.
#' @details The data inputs interactively requested from the user help obtain the following information
#' @details \code{chromosomal segmentation} helps in determining the number of clone/sub-clone cloud to be expected in the data. As
#' variant alleles from different aberrant chromosomes may have similar relative frequencies but discordant clonal
#' interpretation. On the contrary convergent clonal alleles may demonstrate divergent frequencies if arisen from dissimilar aneuploidy.
#' @details \code{clouds} give the program a visual feedback from the user that assume to carry some biological interpretation of the
#' frequency distributions present in the data. This is a subjective estimate that the program later uses for cluster assignment.
#' @details Out of the two methods used for cluster optimization, \emph{GMM} stands for \emph{Gaussian Mixed Models} whereas \emph{bootstrap},
#' as the name suggests perform \emph{bootstrap} resampling of the VAFs in 50 repetitions with 20 runs each to find the most stable parameter
#' for clustering. \emph{GMM} outputs the optimization curve with \code{BIC} or \emph{Bayesian Information Criterion} against number of
#' clusters chosen in the \code{X-axis} where \emph{bootstrap} shows the \code{Smin} statistics instead in the \code{Y-axis}. In both cases the
#' statistics are to be interpreted as proxies for the \emph{entropy} of the system. The maximum entropy is likely to indicate the most stable
#' solution.
#' @details \code{clustering.method} gives the user three choices:
#' @details \code{HKM} is \emph{Heierarchical K-means clustering} which uses heierarchical clustering first to determine the cluster centers
#' that are subsequently used as the starting point for the K-means clustering.
#' \code{bootkm} performs a \emph{bootstrap} resampling of 20 fitted K-means clusters with 50 resamplings to out put the clusters.
#' \code{hybrid} performs \emph{hkm} on the principal component of the data.
#' @details \code{clonality} provides two choices for clonality assignment. The default is \emph{Allelic composition} that measures expected
#' clonality patterns according to the copy numbers. But in cases of unreliable allelic composition estimates this method may fail. In such
#' situations the clonality can be assigned without apriori assumptions with the alternate \emph{density} based method.
#'
#' @return A list of 12 objects is returned that includes all the summary statistics, diagnositics and the predictions as well as the mapping
#' internally used for clonal deconvolution.
#' @return \code{predicted.data} is necessarily an extension to the input \code{data} with the addition of the predicted clone and sub-clone
#' status of each variant for corresponding samples.
#' @return \code{density.map} is a distance matrix convoluted from cluster distances and desity departures.
#' @return \code{collapse} are clusters that are initially prredicted but later collapsed on each other dues to similarity between them.
#' @return \code{fitted.hkm, fitted bootkm or fitted.hybrid} is a vector of initial cluster assignment by the algorithm chosen.
#' Only one of these will have an output and the rest will show \code{NA}.
#' @return \code{Number of unscaled clusters} gives umber of predicted clusters before collapsing with density estimates.
#' @return \code{Number of scaled clusters} gives number of predicted clusters after collapsing (if any).
#' @return \code{cluster.diagnostics} if the optimization method was chosen to be \emph{GMM}, this is an object of \code{S3} class
#' that includes clustering diagnostics from the model-based clustering. If the chosen method was \emph{bootstrap} then this is a list.
#' @return \code{cluster centers} are the centroids of the predicted scaled clusters.
#' @return \code{cluster mapping} provides the map between scaled clusters and the clonal deconvolution assignments
#' @return \code{Dunn index} is the Dunn index for the fitted cluster.
#'
#' @seealso \code{\link{seqn.scale}} \code{\link{cluster.doubt}}
#'
#' @examples
#' \donttest{#cluster.doc(test.dat, 1, 2, optimization.method = 'GMM', clustering.method = 'HKM')}
#' @export
cluster.doc <- function (data = NULL, sample = NULL, vaf = NULL, allele.comp = NULL, n.clone = NULL, n.subclone = NULL, optimization.method = "GMM",
                         clustering.method = "HKM", clonality = "Allelic composition", instruct = TRUE){
  if (missing(sample) | missing(vaf)) stop("Missing entries for sample or vaf")
  if (typeof (sample) != typeof(vaf)) stop("input for sample and vaf should be similar object")
  if (typeof(sample) == "character")
  {x<-as_tibble(data[,c(sample,vaf)])} else {
    x<-as_tibble(subset(data,select=c(sample,vaf)))}
  user_col = colnames(x)
  colnames(x)<-c("sample","vaf")
  vaf.name<-colnames(data)[vaf]

  if(missing(allele.comp) | (missing(n.clone) | missing(n.subclone))){
    print("Insufficient entries provided!")
    i.1 <- Init.input(x, al.comp = allele.comp)
  } else {
    h<-allele.comp
    h1<-unlist(strsplit(h,"+")[1])[1]
    h2<-unlist(strsplit(h,"+")[1])[3]
    h1 <- as.numeric(h1)
    h2 <- as.numeric(h2)
    i.1 <- list(h1=h1, h2=h2, r=as.numeric(n.clone), m=as.numeric(n.subclone))
    g<-dot_plot(d = x)
    g<-g + theme(legend.position="none") + scale_color_brewer(palette="Pastel1")
    print(g)
  }
  tmp.0<-i.1$r+i.1$m


  ######################################################################################################
  ################################# Number of cluster optimization #####################################
  ######################################################################################################

  if (optimization.method == 'GMM') {
    print("Optimizing number of clusters with GMM")
    d_clust <- Mclust(as.matrix(x[,2]), G=2:10)##Deterministic clustering with GMM based density estimate
    m.best <- dim(d_clust$z)[2]                       ##Optimum number of clusters
    cat("Predicted optimal number of clusters:", m.best, "\n")
    print("Plotting BIC curve")
    suppressWarnings(print(plot.bic(d_clust$BIC)))
  } else if (optimization.method == 'bootstrap') {
    print("Optimizing number of clusters with Bootstrapping")
    cat("Bootstrapping is like a box of chocolates. You never know what you're gonna get \n")
    boot.k<-k.select(x[,2], range = 2:10, B = 50, r = 20, threshold = 0.75, scheme_2 = TRUE)
    m.best <- boot.k$k
    cat("Predicted optimal number of clusters:", m.best, "\n")
    print("Plotting Smin curve")
    suppressWarnings(print(plot.bic(boot.k)))
  } else stop("Is that a real optimization method ? .... I felt a great disturbance in the source")

  ######################################################################################################
  ############################ Getting User input about the suggested prediction #######################
  ######################################################################################################

  if (instruct == TRUE){
    if (m.best == tmp.0) {
      print("Your guess was correct! Proceeding with analysis.")
      clust.num = m.best
    } else {
      h <- (readline("Would you like to see my suggestion instead? (Y/N): \n"))
      if (h == "y" | h == "Y") {
        clust.num = m.best
        if((m.best %% 2) == 0) {
          q = m.best/2
          cat(paste0("\n I suggest ", q, " clone(s) and ", q, " subclone(s) \n"))
        } else {
          q = m.best-i.1$t.0
          cat(paste0("\n I suggest ", i.1$t.0, " clone(s) and ", q, " subclone(s) \n"))
        }
      } else if (h == "n" | h == "N"){
        cat("Suspected and predicted clone/subclone segregation mismatch!. \n")
        cat("Proceeding with analysis. \n")
        clust.num = tmp.0
      } else {stop('Fool of a took! You shall not pass!')}
    }
  }  else if (instruct == FALSE){
    cat("Suspected and predicted clone/subclone segregation mismatch!. \n
        Proceeding with analysis. \n")
    clust.num = tmp.0
  } else {stop('Unrecognized input provided!')}
  if (clust.num > 5) {cat("WARNING: More than 5 clusters indicated; may result in redundancy \n")}

  ######################################################################################################
  #################################### Prediction of clusters ##########################################
  ######################################################################################################

  if (clustering.method == 'HKM'){                    ##Predicting Clone Sub-clone status with HKM
    print("Performing Heierarchical K-means clustering")
    hkm<-hkmeans(x[,2],clust.num)
    pred.clust<-hkm$cluster
  } else if (clustering.method == 'bootkm'){            ##Predicting Clone Sub-clone status with bootcluster
    print("Performing Bootstrapped K-means clustering")
    bootkm<-stability(x[,2], clust.num, B=50, r=20,  scheme_2 = TRUE)
    pred.clust<-as.vector(bootkm$membership)
  } else if (clustering.method == 'hybrid'){            ##Predicting Clone Sub-clone status with HCPC
    print("Performing hybrid clustering")
    hybrid<-suppressWarnings(HCPC(x[,2], iter.max=50, nb.clust = clust.num, graph=F))
    pred.clust<-as.numeric(hybrid$data.clust$clust)
  }  else stop("That clustering method... What we've got here is failure to communicate.")

  ######################################################################################################
  ############################## Density based dimension adjustment ####################################
  ######################################################################################################

  tmp<-x
  tmp$clust<-as.numeric(pred.clust)
  tmp.col<-ncol(tmp)
  for (i in 1 :9) { tmp<-cbind(tmp,as.numeric(dbscan(dist(tmp$vaf),0.1+0.02*i)$cluster)) }
  tmp$net<-rowSums(tmp[,(tmp.col+1) : (tmp.col+9)])
  tmp<-tmp[,c(1,2,3,ncol(tmp))]
  tmp1<-as.matrix(table(tmp$clust,tmp$net))

  for(i in 1 : nrow(tmp1)) {
    l<-sum(as.numeric(tmp1[i,]))
    tmp1[i,]<-round(as.numeric(tmp1[i,])/l,2)
  }
  d<-round(dist(tmp1),2)

  a<-vector()
  for(i in 1 : nrow(tmp1)) {a<-c(a,mean(subset(tmp,tmp$clust==i)$vaf))}
  mean.d<-dist(a)
  net.d<-as.matrix(round(d * mean.d,3))               ##Density based scaling of cluster assignments
  diag(net.d)<-1
  similar.d<-sort(as.numeric(rownames(which(net.d< 0.05, arr.ind = T))))
  if (is_empty(similar.d) == FALSE){
    if(length(similar.d)>2) {
      stop("More than one clusters are showing redundancy. I recommend lowering number of clones/sub-clones")}
    else{
      print("Reducing cluster dimension based on density map")
      pred.clust<-replace(pred.clust,pred.clust==similar.d[2],similar.d[1])
    }
  }

  tryCatch({
    if(dim(table(pred.clust))<2) {
      stop("The clusters merged to form only one cloud hence all are predicted clonal")
    }

    df<-scale(x[,2])                                    ##Dunn index for suggested clustering
    clust.stats <- suppressWarnings(cluster.stats(d = dist(df), pred.clust))
    Dunn.index <- as.numeric(clust.stats$dunn)

    ######################################################################################################
    #################################### Clonality assignment ############################################
    ######################################################################################################
    if (clonality == 'density'){
      print("Performing density based clonal deconvolution")
      if (is_empty(similar.d) == FALSE){
        clust.d<-net.d[-similar.d[2],-similar.d[2]]
        a<-a[-similar.d[2]]
      } else {
        clust.d<-net.d
      }

      if (length(a) < 2) {warning("All predicted clusters merged to form only one cloud.
      Consider removing outliers or increase number of clonal / sub-clonal clouds")}

      if(min(dist(a))<0.05){
        warning("CAUTION: centroids of clonal and sub-clonal clusters are too close to comfort (+/- 0.05)!")
      }

      if(min(clust.stats$separation)<0.01){
        warning("CAUTION: minimum separation between points of two different clusters is less than 0.01!")
      }

      clust.a<-as.data.frame(a)
      clust.a$mcgfn<-1:nrow(clust.a)
      rownames(clust.a)<-row.names(clust.d)

      n.cl<-n.sbcl<-0
      sim.pred<-vector()
      cluster.map<-vector()

      for (i in c(1,3,5,7,9)){
        if(nrow(clust.d) > 1){
          next.similar.1<-as.numeric(rownames(which(clust.d == min(clust.d), arr.ind = T)))
          sim.pred[i]<-ifelse(clust.a[which(rownames(clust.a) == next.similar.1[2]),1] >
                                clust.a[which(rownames(clust.a) == next.similar.1[1]),1],'Sub-clone','Clone')
          sim.pred[i+1]<-ifelse(clust.a[which(rownames(clust.a) == next.similar.1[2]),1] <
                                  clust.a[which(rownames(clust.a) == next.similar.1[1]),1],'Sub-clone','Clone')
          cluster.map<-c(cluster.map,next.similar.1)
          n.cl<-n.cl+1
          n.sbcl<-n.sbcl+1
          clust.d<-as.matrix(clust.d[-which(rownames(clust.a) %in% next.similar.1), -which(rownames(clust.a) %in% next.similar.1)])
          clust.a<-clust.a[!rownames(clust.a) %in% next.similar.1,]
        }
        if(nrow(clust.d) < 2) break
      }

      if (nrow(clust.d) == 1) {
        if (i.1$r < i.1$m) {
          sim.pred <- c(sim.pred,'Sub-clone')
        } else {
          sim.pred <- c(sim.pred,'Clone')
        }
        cluster.map <- c(cluster.map, as.numeric(rownames(clust.a)))
      }

      if(i.1$m == 0 | i.1$r == 0) {
        if (i.1$m == 0) {
          sim.pred<-rep("Clone",length(sim.pred))
          cat("Fly you fools! User Suggests there's only Clones \n")
        }
        if (i.1$r == 0) {
          sim.pred<-rep("Sub-clone",length(sim.pred))
          cat("Fly you fools! User Suggests there's only Sub-clones \n")
        }
      }
      pred.map<-as.data.frame(cluster.map)
      pred.map$pred<-sim.pred
      x$cluster.map<-pred.clust
      x<-inner_join(x,pred.map,by='cluster.map')
      x<-x[,c('sample','vaf','pred')]
    }  else {
      print("Allelic composition based clonal deconvolution")

      x$cluster.map<-pred.clust
      mean.vaf <- x %>% group_by(cluster.map) %>% summarize(mean.vaf = round(mean(vaf),3))
      map<-clonality(h1 = i.1$h1, h2 = i.1$h2, center=mean.vaf$mean.vaf)
      map$cluster.map <- mean.vaf[order(mean.vaf[,2], decreasing = TRUE),]$cluster.map
      pred.map<-map[,c(4,3)]
      x<-inner_join(x,pred.map,by='cluster.map')
      x<-x[,c('sample','vaf','clonality')]
      colnames(x)<-c('sample','vaf','pred')
    }

    clone.perc<-table(x$pred)[[1]]*100/nrow(x)
    if ( length(table(x$pred)) > 1 & (clone.perc < 5 | clone.perc > 95)) {
      warning("Too unstable clustering as less than 5% variants are represented in a clone/subclone.
            Recommended: removal of outliers or increase the number of clonal/sub-clonal cloud")
    }

    ######################################################################################################
    #################################### Creating the output #############################################
    ######################################################################################################

    if (vaf.name == 'scaled.vaf') {
      plot.data <- cbind(data[,c(sample,vaf-1)],x$pred)
      colnames(plot.data)<-colnames(x) } else {
        plot.data <- x
      }
    g<-dot_plot(plot.data)
    g<-g + scale_color_manual(values=c("tan1","mediumpurple3"))
    print(g)

    data$'Predicted clonal structure'<-x$pred

    if (exists("d_clust")==TRUE) {cd = d_clust}
    else if (exists("boot.k")==TRUE) {cd = boot.k}
    else {cd = NA}

    invisible(list(predicted.data=data,
                   density.map=net.d,
                   collapse=similar.d,
                   fitted.hkm=ifelse(exists("hkm")==TRUE,hkm,NA),
                   fitted.bootkm=ifelse(exists("bootkm")==TRUE,bootkm,NA),
                   fitted.hybrid=ifelse(exists("hybrid")==TRUE,hybrid,NA),
                   'Number of unscaled clusters' = clust.num,
                   'Number of scaled clusters' = nrow(pred.map),
                   'cluster diagnostics' = cd,
                   'cluster centers' = a,
                   'cluster mapping' = pred.map,
                   'Dunn index' = Dunn.index,
                   'user input' = user_col))
  }, error=function(e){
    cat("ERROR : ",conditionMessage(e), "\n")
    data1=data
    data1$'Predicted clonal structure'<-"clone"
    invisible(list(predicted.data=data1,
                   fitted.hkm=ifelse(exists("hkm")==TRUE,hkm,NA),
                   fitted.bootkm=ifelse(exists("bootkm")==TRUE,bootkm,NA),
                   fitted.hybrid=ifelse(exists("hybrid")==TRUE,hybrid,NA),
                   'Number of unscaled clusters' = clust.num,
                   'user input' = user_col))
  })
}
