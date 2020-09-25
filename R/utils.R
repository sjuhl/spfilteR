
# calculate p-values
pfunc <- function(z,alternative){
  if(alternative=="greater"){
    p <- pnorm(z, lower.tail=F)
  } else if(alternative=="lower"){
    p <- pnorm(z, lower.tail=T)
  } else p <- 2*pnorm(abs(z), lower.tail=F)
  return(p)
}


# calculate empirical p-values
emp.pfunc <- function(draws,z,alternative){
  z <- as.numeric(z)
  # see e.g., North/ Curtis/ Sham (2002) [Am J Hum Genet] for the '+1'
  if(alternative=="greater"){
    p <- (sum(draws>=z)+1)/(length(draws)+1)
  } else if(alternative=="lower"){
    p <- (sum(draws<=z)+1)/(length(draws)+1)
  } else {
    p <- (sum(abs(draws)>=abs(z))+1)/(length(draws)+1) # see e.g., Hartwig (2013) [J Clin Trials]
  }
  return(p)
}


# visual illustration of significance levels
star <- function(p){
  out <- NULL
  out[p<=.001] <- "***"
  out[p<=.01 & p>.001] <- "**"
  out[p<=.05 & p>.01] <- "*"
  out[p<=.1 & p>.05] <- "."
  out[p>.1] <- " "
  return(out)
}


# determine ideal candidate set size (e.g., Chun et al. 2016)
candsetsize <- function(npos, zMI){
  denominator <- 1+exp(2.1480-(6.1808*(zMI+.6)^.1742)/npos^.1298 + 3.3534/(zMI+.6)^.1742)
  nc <- npos/denominator
  return(round(nc,0))
}
