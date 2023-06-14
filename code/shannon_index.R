
# Jialin He edited in 202306


shannon_index <- function(microdat,metadata,group,sample_in_row,
                          p.adj=T,p.signif=T,title=NULL,color=NULL){
  # library(tidyverse)
  # library(ggpubr)
  # library(vegan)
  # library(ggthemes)
  
  if(!require(pacman))install.packages(pacman)
  pacman::p_load(tidyverse,ggpubr,rstatix,vegan,ggthemes,scico)
  
  if(sample_in_row == F){
    microdat <- t(microdat) %>% as.data.frame()
  }
  
  id <- row.names(metadata)
  microdat <- microdat[id,]
  
  ## shannon index ----
  shannon <- data.frame(vegan::diversity(microdat,"shannon"))
  names(shannon) <- "shannon"
  shannon$ID <- row.names(shannon)
  mer <- merge(shannon,metadata, by = 'row.names')
  colnames(mer)[colnames(mer)==group] <- 'group_in_function'
  
  len <- length(unique(mer$group_in_function))
  max <- max(shannon$shannon)
  min <- min(shannon$shannon)
  
  if(is.null(color)){
    color <- scico::scico(n = len,palette = 'berlin')
  }

  label <- ifelse(p.adj,'p.adj','p')
  signif <- ifelse(p.adj,'p.adj.signif','p.signif')
    
  if(len<2){
    stop('ERROR: only one level in your group')
  }else if(len==2){
    stat.test <- compare_means(shannon~group_in_function,
                               p.adjust.method = 'BH',
                               data = mer,method = 'wilcox.test') %>% add_significance()
    if(p.signif){
      label <- signif
    }else{
      stat.test[[label]] <- paste0('P = ',round(stat.test[[label]],3))
    }
    
    shannon_plot <- 
      ggplot(data = mer,aes(x = group_in_function,
                            y = shannon,
                            color = group_in_function))+
      geom_boxplot(width = 0.35,size = 0.66,outlier.shape = NA,
                   alpha = 0.9)+
      geom_jitter(aes(fill = group_in_function),
                  width = 0.25,shape = 21,
                  size = 1,alpha = 0.9)+
      # coord_cartesian(ylim = ylim1*1.05)+
      theme_few()+
      scale_color_manual(values = color)+
      scale_fill_manual(values = color)+
      scale_y_continuous()+
      labs(title = title,x = NULL,y = "Shannon index")+
      theme(plot.title=element_text(family = "serif",
                                    size=12,hjust=0),
            legend.position = 'none',
            axis.title.y= element_text(family = "serif",size=12),
            axis.text.x= element_text(family = "serif",size=12),
            axis.text.y= element_text(family = "serif",size=12),
            axis.ticks.length.x = unit(0.1,'cm'),
            axis.ticks.length.y = unit(0.1,'cm'),
            plot.margin = unit(c(0.2,0.2,0.2,0.2),"cm"))+
      # guides(color=guide_legend(override.aes = list(size=2)))+
      stat_pvalue_manual(data = stat.test,
                         y.position = max+(max-min)/10,
                         family = 'serif',
                         step.increase = 0.09,
                         tip.length = 0.02,label = 'p.adj.signif')
  }else{
    stat.test <- compare_means(shannon~group_in_function,
                               method = 'wilcox.test',
                               data = mer,
                               p.adjust.method = 'BH') %>% add_significance()
    
    if(p.signif){
      label <- signif
    }else{
      stat.test[[label]] <- paste0(rep('P = ',factorial(len)/(2*factorial((len-2)))),
                                   round(stat.test[[label]],3))
    }
    
    shannon_plot <- 
      ggplot(data = mer,aes(x = group_in_function,
                            y = shannon,
                            color = group_in_function))+
      geom_boxplot(width = 0.35,size = 0.66,outlier.shape = NA,
                   alpha = 0.9)+
      geom_jitter(aes(fill = group_in_function),
                  width = 0.25,shape = 21,
                  size = 1,alpha = 0.9)+
      # coord_cartesian(ylim = ylim1*1.05)+
      theme_few()+
      scale_color_manual(values = color)+
      scale_fill_manual(values = color)+
      scale_y_continuous()+
      labs(title = title,x = NULL,y = "Shannon index")+
      theme(plot.title=element_text(family = "serif",
                                    size=12,hjust=0),
            legend.position = 'none',
            axis.title.y= element_text(family = "serif",size=12),
            axis.text.x= element_text(family = "serif",size=12),
            axis.text.y= element_text(family = "serif",size=12),
            axis.ticks.length.x = unit(0.1,'cm'),
            axis.ticks.length.y = unit(0.1,'cm'),
            plot.margin = unit(c(0.2,0.2,0.2,0.2),"cm"))+
      stat_pvalue_manual(data = stat.test,
                         y.position = max+(max-min)/10,
                         step.increase = 0.09,
                         tip.length = 0.02,
                         label = label,
                         family = 'serif',size = 3.8)
  }
  list <- list(data=shannon,plot=shannon_plot)
  return(list)
}
