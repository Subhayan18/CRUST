pqn <- function (sample, ref = NULL){
  S.sample = rowSums(abs(sample), na.rm = T)
  sample = sample/S.sample
  N.ref = apply(ref, 2, median, na.rm = T)
  N.sample = t(t(sample)/N.ref)
  param = apply(N.sample, 1, median, na.rm = T)
  return(N.sample = sample/param)
}
