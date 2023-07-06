# Jialin He 202307

mantel_function <- function(microdat,
                            metadata,
                            sample_in_row=T,
                            var,
                            method = 'pearson',
                            community_name='Fungi',
                            permutations = 999,
                            parallel = 4,
                            seed = 0,
                            rbreak = c(-Inf,0.1,0.2,Inf),
                            rlabel = c("<0.1","0.1-0.2",">0.2"),
                            pbreak = c(-Inf,0.001,0.01,0.05,Inf),
                            plabel = c("<0.001","<0.01","<0.05",">0.05")){
  if(!require(pacman))install.packages(pacman)
  pacman::p_load(tidyverse,vegan)
  
  if(sample_in_row == F){
    microdat <- t(microdat) %>% as.data.frame()
  }
  
  metadata <- metadata[row.names(metadata) %in% row.names(microdat),]
  id <- row.names(metadata)
  microdat <- microdat[id,]
  
  r <- NULL
  p <- NULL
  for(i in var){
    meta <- metadata[!is.na(metadata[[i]]),]
    micro <- microdat[row.names(meta),]
    
    dist1 <- vegdist(micro,method = 'bray')
    dist2 <- vegdist(meta[[i]],method = 'euclidean')
    
    set.seed(seed)
    man <- mantel(dist1,dist2,method = method,
                  permutations = permutations,
                  parallel = parallel)
    
    r <- c(r,man$statistic)
    p <- c(p,man$signif)
    
  }
  
  dat <- data.frame(name=rep(community_name,length(var)),
                    variable=var,
                    r=r,
                    p=p)
  # 以我的经验,r<0的时候就不存在P<0.05的情况
  dat$rlabel <- cut(dat$r,breaks = rbreak,
                    labels = rlabel,
                    include.lowest = TRUE,right = FALSE)
  
  # 主要是因为permutation 999，没办法得到<0.001的P
  # 所以right = T
  dat$plabel <- cut(dat$p,breaks = pbreak,
                    labels = plabel,
                    include.lowest = F,right = T)
  return(dat)
}

