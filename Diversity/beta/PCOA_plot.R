# Jialin He edited in 202306

PCOA_plot <- function(microdat,metadata,group,
                      sample_in_row = T,
                      signif_method = 'adonis',
                      distance = 'bray',
                      parallel =4,
                      title = NULL,
                      color = NULL,
                      legend.position = 'top',
                      path = NULL,
                      filename = 'PCoA',
                      seed = 0,
                      width = 6.6,
                      height = 5.1){
  if(!require(pacman))install.packages(pacman)
  pacman::p_load(tidyverse,ggpubr,vegan,car,ggthemes,scico)

  if(sample_in_row == F){
    microdat <- t(microdat) %>% as.data.frame()
  }
  
  metadata <- metadata[row.names(metadata) %in% row.names(microdat),]
  id <- row.names(metadata)
  microdat <- microdat[id,]

  colnames(metadata)[colnames(metadata)==group] <- 'group_in_function'
  
  len <- length(unique(metadata$group_in_function))
  
  if(is.null(color)){
    color <- scico::scico(n = len,palette = 'berlin')
  }
  
  if(len<2){
    stop('ERROR: only one level in your group')
  }else{
    # calculate P value
    
    set.seed(seed)
    if(signif_method == 'anosim'){
      
      anosim <- anosim(microdat,metadata$group_in_function,
                       permutations = 999,distance = distance,
                       parallel = parallel)
      P_val <- anosim$signif %>% round(.,3)
      statistic <- anosim$statistic %>% round(.,3)
      label <- bquote(atop(R ==.(statistic),
                           italic(P) ==.(P_val)))
    }else if(signif_method == 'adonis'){
      
      adonis <- adonis2(microdat~group_in_function,data = metadata,
                        permutations = 999,method = distance,
                        parallel = parallel)
      P_val <- adonis[1,5] %>% round(.,3)
      statistic <- (adonis[1,3]*100) %>% round(.,3)
      statistic <- paste0(statistic,'%')
      label <- bquote(atop(R^2 ==.(statistic),
                           italic(P) ==.(P_val)))
    }
    
    # PCOA
    # distance matrix
    pcoa_dis <- vegdist(microdat, method = distance)
    
    pcoa <- cmdscale(pcoa_dis, k = (nrow(microdat) - 1), eig = TRUE)
    
    pcoa_exp <- pcoa$eig/sum(pcoa$eig)
    pcoa1 <- paste('PCoA1 :', round(100*pcoa_exp[1], 2), '%')
    pcoa2 <- paste('PCoA2 :', round(100*pcoa_exp[2], 2), '%')
    
    PC <- data.frame(pcoa$point)[1:2]
    PC <- cbind.data.frame(PC,metadata$group_in_function)
    names(PC) <- c('pcoa1','pcoa2','group_in_function')
    
    group_average <- aggregate(cbind(pcoa1, pcoa2)~group_in_function, 
                               data = PC, FUN = mean)
    
    xmin <- NULL
    xmax <- NULL
    ymin <- NULL
    ymax <- NULL
    
    for(i in unique(metadata$group_in_function)){
      temp <- PC[PC$group_in_function == i,]
      eli <- dataEllipse(x = temp$pcoa1,y = temp$pcoa2,
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
  }
  
  plot <- ggplot(PC)+
    geom_point(aes(x=pcoa1, y=pcoa2,
                   color = group_in_function),
               size = 0.8)+
    stat_ellipse(aes(x = pcoa1, y = pcoa2, 
                     color = group_in_function),
                 level = 0.95, linetype = 2,linewidth = 0.8,
                 show.legend = FALSE) +
    geom_point(data = group_average, 
               aes(x = pcoa1, y = pcoa2, 
                   color = group_in_function),
               size = 2.2, show.legend = FALSE)+
    scale_x_continuous(limits = c(xmin-(xmax-xmin)*0.025,
                                  xmax+(xmax-xmin)*0.025))+
    scale_y_continuous(limits = c(ymin-(ymax-ymin)*0.025,
                                  ymax+(ymax-ymin)*0.025))+
    guides(color=guide_legend(override.aes = list(size=2.5)))+
    labs(title = title,x = pcoa1,y = pcoa2,color = '')+
    scale_color_manual(values = color)+
    annotate('text', label = label, 
             x = xmax*0.9, y = ymax*0.95, 
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

  names(PC) <- c('PCoA1','PCoA2',group)
  
  list <- list(plot = plot,
               dat = PC)
  
  return(list)
}


