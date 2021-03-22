seg.track <- function(a,b){
  b$Bf<-as.numeric(b$Bf)
  b$position <- as.numeric(b$position)
  b$depth.ratio <-as.numeric(b$depth.ratio)
  AI <- vector()
  R <- vector()
  for (i in 1 : nrow(a)){
    test<-subset(b, b$chromosome == a$chromosome[i] & (b$position >= a$start.pos[i] & b$position <= a$end.pos[i]))
    hkm<-factoextra::hkmeans(abs(test$Bf-0.5),2)
    AI <- c(AI,min(hkm$centers) / max(hkm$centers))
    R <- c(R,mean(log(subset(b, b$chromosome == a$chromosome[i] &
                               (b$position >= a$start.pos[i] & b$position <= a$end.pos[i]))$depth.ratio,2)))
  }
  list(AI,R)
}
