
# Jialin He edited in 202306


alpha_diversity <- function(microdat,
                            metadata,group,
                            sample_in_row,
                            plot_index = 'Shannon',
                            p.adj = T,
                            p.signif = T,
                            title = NULL,
                            color = NULL,
                            path = NULL,
                            filename = 'diversity',
                            width = 6.6,
                            height = 5){
  
  if(!require(pacman))install.packages(pacman)
  pacman::p_load(tidyverse,ggpubr,rstatix,vegan,scico)
  
  if(sample_in_row == F){
    microdat <- t(microdat) %>% as.data.frame()
  }
  
  # 去除所有样本都为0的OTU/ASV/Species/Genus etc.
  microdat <- microdat[,colSums(microdat)>0]
  
  id <- row.names(metadata)
  microdat <- microdat[id,]
  
  ## α-diversity ----
  Shannon <- vegan::diversity(microdat,"shannon")
  Simpson <- vegan::diversity(microdat,'simpson')
  Gini_Simpson <- 1-vegan::diversity(microdat,index = 'simpson')
  Inv_Simpson <- vegan::diversity(microdat,index = 'invsimpson')
  Richness <- apply(microdat,1,function(x){sum(x>0)})
  # Richness <- specnumber(microdat)
  
  # evenness
  Pielou <- Shannon/log(Richness,base = exp(1))
  equitability <- Inv_Simpson/Richness
  
  # AVD
  ai <- abs(microdat-colMeans(microdat))/apply(microdat,2,sd)
  AVD <- rowSums(ai)/(1*ncol(microdat))
  
  index <- cbind(Shannon,Simpson,Gini_Simpson,Inv_Simpson,Richness,
                 Pielou,equitability,AVD) %>% as.data.frame()
  
  mer <- merge(index,metadata, by = 'row.names') %>% 
    column_to_rownames(.,'Row.names')
  colnames(mer)[colnames(mer)==group] <- 'group_in_function'
  mer <- mer[,c('Shannon','Simpson','Gini_Simpson','Inv_Simpson','Richness',
                'Pielou','equitability','AVD','group_in_function')]
  index <- mer %>% rename(.,group=group_in_function)
  
  if(plot_index == 'Shannon'){
    ytitle <- 'Shannon Index'
  }else if(plot_index == 'Gini_Simpson'){
    ytitle <- 'Gini-Simpson Index'
  }else if(plot_index == 'Simpson'){
    ytitle <- 'Simpson Index'
  }else if(plot_index == 'Inv_Simpson'){
    ytitle <- 'Inverse Simpson Index'
  }else if(plot_index == 'Richness'){
    ytitle <- 'Richness'
  }else if(plot_index == 'Pielou'){
    ytitle <- 'Pielou`s Evenness'
  }else if(plot_index == 'equitability'){
    ytitle <- 'Simpson`s Evenness'
  }else if(plot_index == 'AVD'){
    ytitle <- 'AVD Index'
  }else{
    stop('ERROR: This index did not exist in it')
  }
  
  len <- length(unique(mer$group_in_function))
  max <- max(mer[[plot_index]])
  min <- min(mer[[plot_index]])
  
  if(is.null(color)){
    color <- scico::scico(n = len,palette = 'berlin')
  }
  
  label <- ifelse(p.adj,'p.adj','p')
  signif <- ifelse(p.adj,'p.adj.signif','p.signif')
  
  if(len<2){
    stop('ERROR: only one level in your group')
  }else{
    colnames(mer)[colnames(mer)==plot_index] <- 'y'
    stat.test <- compare_means(y~group_in_function,
                               p.adjust.method = 'BH',
                               data = mer,method = 'wilcox.test') %>% 
      add_significance()
    
    if(p.signif){
      label <- signif
    }else{
      stat.test[[label]] <- paste0('P = ',round(stat.test[[label]],3))
    }
    
    plot <- 
      ggplot(data = mer,aes(x = group_in_function,
                            y = y,
                            color = group_in_function))+
      geom_boxplot(width = 0.35,size = 0.66,outlier.shape = NA,
                   alpha = 0.9)+
      geom_jitter(aes(fill = group_in_function),
                  width = 0.25,shape = 21,
                  size = 1,alpha = 0.9)+
      # coord_cartesian(ylim = ylim1*1.05)+
      theme_bw()+
      scale_color_manual(values = color)+
      scale_fill_manual(values = color)+
      scale_y_continuous()+
      labs(title = title,x = NULL,y = ytitle)+
      theme(plot.title=element_text(family = "serif",
                                    size=15,hjust=0),
            legend.position = 'none',
            axis.title.y= element_text(family = "serif",size=15),
            axis.text.x= element_text(family = "serif",size=12),
            axis.text.y= element_text(family = "serif",size=12),
            axis.ticks.length.x = unit(0.1,'cm'),
            axis.ticks.length.y = unit(0.1,'cm'),
            plot.margin = unit(c(0.2,0.2,0.2,0.2),"cm"))+
      # guides(color=guide_legend(override.aes = list(size=2)))+
      stat_pvalue_manual(data = stat.test,
                         y.position = max+(max-min)/10,
                         family = 'serif',
                         step.increase = 0.09,size = 3.8,
                         tip.length = 0.02,label = label)
  }
  
  list <- list(data=index,plot=plot)
  
  if(!is.null(path)){
    ggsave(plot = plot,
           filename = paste0(path,'/',filename,'_plot.pdf'),
           width = width,height = height)
    index <- rownames_to_column(index,'ID')
    write.table(index,paste0(path,'/alpha_diversity.txt'),
                quote = F,row.names = F,sep = '\t')
  }
  
  return(list)
}
