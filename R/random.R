#' Function to easier call runif
crunif <- function(invec) {
  
  if (length(invec) != 2L) {
    stop("wrong size")
  }
  
  runif(1L, invec[1], invec[2])
}

#' Function to easier call rnorm
crnorm <- function(invec) {
  
  if (length(invec) != 2L) {
    stop("wrong size")
  }
  
  # 95% CI:
  # CI: m +/- q*sd
  # q <- 1.959964 # 95%
  # invec[1] = m - q*sd
  # sd = (m - invec[1])/q
  # 
  m <- mean(invec)
  s <- (m - invec[1])/1.959964
  
  if (FALSE) {
    qnorm(0.025, m, s)
    qnorm(0.975, m, s)
  }
  
  x <- rnorm(1L, m, s)
  
  while (x < 0) {
    x <- rnorm(1L, m, s)
  }
  
  # if (x < 0) {
  #   print(invec)
  #   print(m)
  #   print(s)
  #   print(qnorm(0.025, m, s))
  #   print(qnorm(0.975, m, s))
  #   print(x)
  #   stop("x < 0 not expected")
  # }
  
  return(x)
}
