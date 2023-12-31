# Jialin He 202306

NMDS <- function(microdat,
                 metadata,
                 group,
                 distance = 'bray',
                 k = 2,try = 20,trymax = 50,
                 autotransform = T,
                 parallel = 4,
                 signif_method = 'adonis',
                 title = NULL,
                 color = NULL,
                 legend.position = 'top',
                 path = NULL,
                 filename = 'NMDS',
                 seed = 0,
                 width = 6.6,
                 height = 5.1){
  
  if(!require(pacman))install.packages(pacman)
  pacman::p_load(tidyverse,ggpubr,vegan,car,ggthemes,scico)

  if(sum(sapply(microdat,is.numeric))!=ncol(microdat)){
    stop('ERROR: Only numeric values can be included in the microdat')
  }
  
  if(sum(row.names(metadata) %in% row.names(microdat))>0){
    sample_in_row = T
    print('Sample ID in row')
  }else if(sum(row.names(metadata) %in% colnames(microdat))>0){
    sample_in_row = F
    print('Sample ID in col')
  }else{
    stop('ERROR: Sample ID must be the rowname/colname in both dataset')
  }
  
  if(sample_in_row == F){
    microdat <- t(microdat) %>% as.data.frame()
  }
  
  metadata <- metadata[row.names(metadata) %in% row.names(microdat),]
  id <- row.names(metadata)
  microdat <- microdat[id,]
  
  colnames(metadata)[colnames(metadata) == group] <- 'group_in_function'
  
  set.seed(seed)
  nmds <- metaMDS(microdat,distance = distance,k = k,
                  try = try,trymax = trymax,
                  autotransform = autotransform,
                  model = 'global',noshare = F,
                  trace = 2,parallel = parallel)
  
  
  # stressplot(nmds)
  # goodness <- goodness(nmds)
  # plot(nmds)
  # points(nmds,cex = goodness*100)
  
  nmds_dat <- as.data.frame(scores(nmds)$site)
  dat <- cbind(metadata,nmds_dat)
  len <- length(unique(dat$group_in_function))
  
  if(len<2){
    stop('ERROR: only one level in your group')
  }
  
  if(is.null(color)){
    color <- scico::scico(n = len,palette = 'berlin')
  }
  
  xmin <- NULL
  xmax <- NULL
  ymin <- NULL
  ymax <- NULL
  
  for(i in unique(dat$group_in_function)){
    temp <- dat[dat$group_in_function == i,]
    eli <- dataEllipse(x = temp$NMDS1,y = temp$NMDS2,
                       levels = 0.95,draw = F) %>% as.data.frame()
    xmint <- min(eli$x)
    xmaxt <- max(eli$x)
    ymint <- min(eli$y)
    ymaxt <- max(eli$y)
    
    xmin <- c(xmin,xmint)
    xmax <- c(xmax,xmaxt)
    ymin <- c(ymin,ymint)
    ymax <- c(ymax,ymaxt)
  }
  
  xmin <- min(xmin)
  xmax <- max(xmax)
  ymin <- min(ymin)
  ymax <- max(ymax)
  
  
  set.seed(seed)
  if(signif_method == 'anosim'){
    
    anosim <- anosim(microdat,dat$group_in_function,
                     permutations = 999,distance = distance,
                     parallel = parallel)
    P_val <- anosim$signif %>% round(.,3)
    statistic <- anosim$statistic %>% round(.,3)
    
  }else if(signif_method == 'adonis'){
    
    adonis <- adonis2(microdat~group_in_function,data = dat,
                      permutations = 999,method = distance,
                      parallel = parallel)
    P_val <- adonis[1,5] %>% round(.,3)
    statistic <- (adonis[1,3]*100) %>% round(.,3)
    statistic <- paste0(statistic,'%')
    
  }
  
  stress <- round(nmds$stress,3)
  label <- bquote(atop(italic(Stress) ==.(stress),
                       italic(P) ==.(P_val)))
  
  plot <- ggplot(dat)+
    geom_point(aes(x=NMDS1, y=NMDS2,
                   color = group_in_function),
               size = 0.8)+
    stat_ellipse(aes(x = NMDS1, y = NMDS2, 
                     color = group_in_function),
                 level = 0.95, linetype = 2,linewidth = 0.8,
                 show.legend = FALSE) +
    scale_x_continuous(limits = c(xmin-(xmax-xmin)*0.025,
                                  xmax+(xmax-xmin)*0.025))+
    scale_y_continuous(limits = c(ymin-(ymax-ymin)*0.025,
                                  ymax+(ymax-ymin)*0.025))+
    guides(color=guide_legend(override.aes = list(size=2.5)))+
    labs(title = title,color = '')+
    scale_color_manual(values = color)+
    annotate('text',label = label,
             x = xmax*0.88, y = ymax*0.92, 
             size = 4.2, family = 'serif')+
    theme_bw()+
    theme(legend.key = element_blank(),
          legend.position = legend.position,
          plot.title = element_text(family = "serif",size=15),
          axis.title.x= element_text(family = "serif",size=15),
          axis.title.y= element_text(family = "serif",size=15),
          axis.text.x= element_text(family = "serif",size=12),
          axis.text.y= element_text(family = "serif",size=12),
          legend.text = element_text(family = "serif",size=10),
          axis.ticks.length.x = unit(0.1,'cm'),
          axis.ticks.length.y = unit(0.1,'cm'),
          plot.margin = unit(c(0.2,0.2,0.2,0.2),"cm"))
  
  if(!is.null(path)){
    ggsave(plot = plot,
           filename = paste0(path,'/',filename,'_plot.pdf'),
           width = width,height = height)
  }
  
  return(plot)
}

