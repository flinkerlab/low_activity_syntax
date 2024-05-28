### Fisher-z transform (for statistical analyses of correlation r-values)
### July 2023
### adam.milton.morgan@gmail.com

## Transforms correlation r-values (-1 to 1) to normally distributed z-values (-inf to inf)
fisher.z.transform <- function(rs) {
  zs <- 0.5 * log((1 + rs) / (1 - rs))
  return(zs)
}

## Transforms correlation r-values (-1 to 1) to normally distributed z-values (-inf to inf)
undo.fisher.z.transform <- function(zs) {
  rs <- (exp(2 * zs) - 1) / (exp(2 * zs) + 1)
  return(rs)
}

