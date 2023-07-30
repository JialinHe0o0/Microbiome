# Jialin He 202307

abun_filter <- function(dat,
                        sample_in_row = T,
                        count = F,
                        abundance = 0.001){
  
  if(!require(pacman))install.packages(pacman)
  pacman::p_load(tidyverse)
  
  if(sample_in_row == F){
    dat <- t(dat) %>% as.data.frame()
  }
  
  if(count == T){
    dat <- (dat/rowSums(dat)) %>% as.data.frame()
  }
  
  dat <- dat[,colMeans(dat)>=abundance]
  
  return(dat)
}

