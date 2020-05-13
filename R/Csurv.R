C.surv <- function(q, Y, base) {
  result <- prod(base[Y <= q])
  return(result)
}