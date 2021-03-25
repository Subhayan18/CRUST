ChromSeg<-function (sequenza.extract, cp.table = NULL,
               cellularity = NULL, ploidy = NULL, CNt.max = 20, XY = c(X = "X", Y = "Y"),
               chromosome.list = 1:24)
{
  seg.tab <- do.call(rbind, sequenza.extract$segments[chromosome.list])
  seg.len <- (seg.tab$end.pos - seg.tab$start.pos)/1000000
  avg.depth.ratio <- sequenza.extract$avg.depth.ratio

    cint <- get.ci(cp.table)
    if (!is.null(cellularity) || !is.null(ploidy)) {
      if (is.null(cellularity)) {
        cellularity <- cint$max.cellularity
      }
      if (is.null(ploidy)) {
        ploidy <- cint$max.ploidy
      }
      points(x = ploidy, y = cellularity, pch = 5)
      text(x = ploidy, y = cellularity, labels = "User selection",
           pos = 3, offset = 0.5)
    }
    else {
      cellularity <- cint$max.cellularity
      ploidy <- cint$max.ploidy
    }
  mut.tab <- na.exclude(do.call(rbind, sequenza.extract$mutations[chromosome.list]))

    segs.is.xy <- seg.tab$chromosome %in% XY
    mut.is.xy <- mut.tab$chromosome %in% XY

  avg.sd.ratio <- sum(seg.tab$sd.ratio * seg.tab$N.ratio, na.rm = TRUE)/sum(seg.tab$N.ratio,
                                                                            na.rm = TRUE)
  avg.sd.Bf <- sum(seg.tab$sd.BAF * seg.tab$N.BAF, na.rm = TRUE)/sum(seg.tab$N.BAF,
                                                                     na.rm = TRUE)
  cn.alleles <- baf.bayes(Bf = seg.tab$Bf[!segs.is.xy], CNt.max = CNt.max,
                          depth.ratio = seg.tab$depth.ratio[!segs.is.xy], cellularity = cellularity,
                          ploidy = ploidy, avg.depth.ratio = avg.depth.ratio, sd.ratio = seg.tab$sd.ratio[!segs.is.xy],
                          weight.ratio = seg.len[!segs.is.xy], sd.Bf = seg.tab$sd.BAF[!segs.is.xy],
                          weight.Bf = 1, ratio.priority = F, CNn = 2)
  seg.res <- cbind(seg.tab[!segs.is.xy, ], cn.alleles)

  return(dplyr::as_tibble(seg.res)[,c(1,2,3,9,10,11,12)])
}
