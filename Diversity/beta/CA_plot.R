# Jialin He 202308

CA_plot <- function(microdat,
                    metadata,
                    group,
                    color=NULL,
                    scaling=1,
                    seed=0,
                    signif_method='adonis',
                    distance='chi',
                    N=10,
                    parallel=4,
                    legend.position='top',
                    title=NULL){
  
  if(!require(pacman))install.packages(pacman)
  pacman::p_load(tidyverse,ggpubr,vegan,car,ggrepel,ggthemes,scico)
  
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
  colnames(metadata)[colnames(metadata)==group] <- 'group_in_function'
  
  len <- length(unique(metadata$group_in_function))
  
  if(is.null(color)){
    color <- scico::scico(n = len,palette = 'berlin')
  }
  
  if(len<2)stop('ERROR: only one level in your group')
  
  ca <- cca(microdat)
  
  ca_eig1 <- (100*ca$CA$eig[1]/sum(ca$CA$eig)) %>%
    round(2)
  ca_eig2 <- (100*ca$CA$eig[2]/sum(ca$CA$eig)) %>%
    round(2)
  
  ca_sample <- scores(ca,choices = 1:length(ca$CA$eig),scaling = scaling)$sites %>%
    as.data.frame() %>% cbind(metadata$group_in_function)
  colnames(ca_sample)[length(ca$CA$eig)+1] <- 'group_in_function'
  
  ca_feature <- scores(ca,choices = 1:length(ca$CA$eig),scaling = scaling)$species %>%
    as.data.frame()
  
  xtitle <- paste('CA1 :', ca_eig1, '%')
  ytitle <- paste('CA2 :', ca_eig2, '%')
  
  xmean <- aggregate(ca_sample$CA1~group_in_function, 
                     data = ca_sample, FUN = mean)
  ymean <- aggregate(ca_sample$CA2~group_in_function, 
                     data = ca_sample, FUN = mean)
  group_average <- left_join(xmean,ymean,"group_in_function")
  colnames(group_average)[2:3] <- c('x','y')
  
  xmin <- NULL
  xmax <- NULL
  ymin <- NULL
  ymax <- NULL
  
  for(i in levels(metadata$group_in_function)){
    temp <- ca_sample[ca_sample$group_in_function == i,]
    eli <- dataEllipse(x = temp$CA1,y = temp$CA2,
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
  
  
  topN <- colnames(microdat)[colSums(microdat) %>% order(decreasing = T)][0:N]
  topN <- ca_feature[topN,] %>% rownames_to_column('taxo')
  
  plot <- ggplot(ca_sample)+
    geom_point(aes(x=CA1, y=CA2,
                   color = group_in_function),
               size = 0.8)+
    geom_point(data = topN,aes(x = CA1,y = CA2),shape = 17,
               size=2.5,alpha = 0.66)+
    stat_ellipse(aes(x = CA1, y = CA2, 
                     color = group_in_function),
                 level = 0.95, linetype = 2,linewidth = 0.8,
                 show.legend = FALSE)+
    geom_point(data = group_average, 
               aes(x = x, y = y, 
                   color = group_in_function),
               size = 2.2, show.legend = FALSE)+
    geom_text_repel(data = topN,aes(x = CA1,y = CA2,
                                    label = taxo),
                    min.segment.length = 1,
                    #vjust = 2, 
                    hjust=0,max.overlaps = Inf,
                    force = 3,
                    force_pull = 2,
                    box.padding = 0.1,
                    segment.linetype = 1,
                    segment.alpha = 0.8,
                    direction = 'y',
                    family = 'serif')+
    scale_x_continuous(limits = c(xmin-(xmax-xmin)*0.025,
                                  xmax+(xmax-xmin)*0.025))+
    scale_y_continuous(limits = c(ymin-(ymax-ymin)*0.025,
                                  ymax+(ymax-ymin)*0.025))+
    guides(color=guide_legend(override.aes = list(size=2.5)))+
    labs(title = title,x = xtitle,y = ytitle,color = '')+
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
  
  return(plot)
  
}

