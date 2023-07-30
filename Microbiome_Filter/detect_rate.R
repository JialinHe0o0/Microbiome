# Jialin He 202307

detect_rate <- function(dat,
                        sample_in_row = T,
                        proportion = 0.05){

  if(!require(pacman))install.packages(pacman)
  pacman::p_load(tidyverse)
  
  if(sample_in_row == F){
    dat <- t(dat) %>% as.data.frame()
  }
  
  a <- function(x){
    detect <- sum(x>0)/length(x)
    return(detect)
  }
  rate <- sapply(dat,a)
  dat <- dat[,rate>=proportion]
  
  list <- list(rate=rate,
               dat=dat)
  return(list)
}

