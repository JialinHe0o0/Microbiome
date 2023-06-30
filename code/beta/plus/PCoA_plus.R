# Jialin He 202306

PCOA_plus <- function(microdat,metadata,group,
                      sample_in_row = T,
                      signif_method = 'adonis',
                      distance = 'bray',
                      parallel =4,
                      color = NULL,
                      legend.position = c(0.9,0.2),
                      p.adj = 'BH',
                      parametric_test = T,
                      path = NULL,
                      filename = 'PCoA_plus',
                      seed = 0,
                      width = 6.6,
                      height = 5.1){
  if(!require(pacman))install.packages(pacman)
  pacman::p_load(tidyverse,ggpubr,vegan,car,ggthemes,scico,
                 agricolae,patchwork)
  
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
    xtitle <- paste('PCoA1 :', round(100*pcoa_exp[1], 2), '%')
    ytitle <- paste('PCoA2 :', round(100*pcoa_exp[2], 2), '%')
    
    PC <- data.frame(pcoa$point)[1:2]
    PC <- cbind.data.frame(PC,metadata$group_in_function)
    names(PC) <- c('PCoA1','PCoA2','group_in_function')
    
    group_average <- aggregate(cbind(PCoA1, PCoA2)~group_in_function, 
                               data = PC, FUN = mean)
    
    # 保证点的置信区域可以被完整展示
    # 分组计算每个组置信区域在x,y轴最大最小值
    # 取组间最大最小值
    
    xmin <- NULL
    xmax <- NULL
    ymin <- NULL
    ymax <- NULL
    
    for(i in unique(metadata$group_in_function)){
      temp <- PC[PC$group_in_function == i,]
      eli <- dataEllipse(x = temp$PCoA1,y = temp$PCoA2,
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
  
  xmin <- xmin-(xmax-xmin)*0.025
  xmax <- xmax+(xmax-xmin)*0.025
  
  ymin <- ymin-(ymax-ymin)*0.025
  ymax <- ymax+(ymax-ymin)*0.025
  
  # 参数or非参数检验
  
  if(parametric_test){
    model <- aov(PCoA1~group_in_function,data = PC)
    test1 = LSD.test(model,trt = 'group_in_function',
                     p.adj = c(p.adj),group = T) %>% 
      .$group %>% as.data.frame() %>% 
      rownames_to_column('Group.1')
    
    model <- aov(PCoA2~group_in_function,data = PC)
    test2 = LSD.test(model,trt = 'group_in_function',
                     p.adj = c(p.adj),group = T) %>% 
      .$group %>% as.data.frame() %>% 
      rownames_to_column('Group.1')
  }else{
    test1 = agricolae::kruskal(y = PC$PCoA1,
                               trt = PC$group_in_function,
                               p.adj = p.adj) %>% 
      .$group %>% as.data.frame() %>% 
      rownames_to_column('Group.1')
    
    test2 = agricolae::kruskal(y = PC$PCoA2,
                               trt = PC$group_in_function,
                               p.adj = p.adj) %>% 
      .$group %>% as.data.frame() %>% 
      rownames_to_column('Group.1')
  }
  
  
  # 取箱式图最大值 (不是max，因为max可能是离群点)
  
  P100 <- function(x){
    P100 <- boxplot.stats(x)$stats[[5]]
    return(P100)
  }
  
  x100 <- aggregate(PC$PCoA1,list(PC$group_in_function),P100) %>% 
    as.data.frame()
  y100 <- aggregate(PC$PCoA2,list(PC$group_in_function),P100) %>%
    as.data.frame()
  
  test1 <- merge(x100,test1,by = 'Group.1')
  test2 <- merge(y100,test2,by = 'Group.1')
  test1$Group.1 <- factor(test1$Group.1,
                          levels = levels(metadata$group_in_function))
  test2$Group.1 <- factor(test2$Group.1,
                          levels = levels(metadata$group_in_function))
  
  # PCoA1箱式图
  
  plot2 <- ggplot()+
    geom_boxplot(data = PC,aes(x = group_in_function,
                               y = PCoA1,
                               color = group_in_function),
                 width = 0.35,size = 0.66,outlier.shape = NA,
                 alpha = 0.9)+
    coord_flip()+
    # coord_cartesian(ylim = ylim1*1.05)+
    theme_bw()+
    scale_color_manual(values = color)+
    scale_fill_manual(values = color)+
    scale_y_continuous(limits = c(xmin,xmax))+
    scale_x_discrete(limits = rev(levels(PC$group_in_function)))+
    labs(x = '',y = '')+
    geom_text(data = test1,aes(x = Group.1,
                               y = x+(xmax-max(test1$x))/1.5,
                               label = groups,
                               color = Group.1),
              show.legend = F)+
    theme(plot.title=element_text(family = "serif",
                                  size=15,hjust=0),
          legend.position = 'none',
          plot.background = element_blank(),
          axis.title.y= element_text(family = "serif",size=15),
          axis.text.x= element_blank(),
          axis.text.y= element_text(family = "serif",size=10),
          axis.ticks.length.x = unit(0,'cm'),
          axis.ticks.length.y = unit(0.1,'cm'),
          plot.margin = unit(c(0,0,0,0),"cm"))
  
  # PCoA2箱式图
  
  plot3 <- ggplot()+
    geom_boxplot(data = PC,aes(x = group_in_function,
                               y = PCoA2,
                               color = group_in_function),
                 width = 0.35,size = 0.66,outlier.shape = NA,
                 alpha = 0.9)+
    theme_bw()+
    scale_color_manual(values = color)+
    scale_fill_manual(values = color)+
    scale_y_continuous(limits = c(ymin,ymax))+
    labs(title = '',x = '',y = '')+
    geom_text(data = test2,aes(x = Group.1,
                               y = x+(ymax-max(test2$x))/1.5,
                               label = groups,
                               color = Group.1),
              show.legend = F)+
    theme(plot.title=element_text(family = "serif",
                                  size=15,hjust=0),
          legend.position = 'none',
          plot.background = element_blank(),
          axis.title.y= element_text(family = "serif",size=15),
          axis.text.x= element_text(family = "serif",size=10,angle = 45,
                                    vjust = 1,hjust = 1),
          axis.text.y= element_blank(),
          axis.ticks.length.x = unit(0.1,'cm'),
          axis.ticks.length.y = unit(0,'cm'),
          plot.margin = unit(c(0,0,0,0),"cm"))
  
  # temp <- ggplot_build(plot = plot2)
  # xbox <- temp$layout$panel_scales_y[[1]]$range$range
  # temp <- ggplot_build(plot = plot3)
  # ybox <- temp$layout$panel_scales_y[[1]]$range$range
  
  # beta diversity
  
  plot1 <- ggplot(PC)+
    geom_point(aes(x=PCoA1, y=PCoA2,
                   color = group_in_function),
               size = 0.8)+
    stat_ellipse(aes(x = PCoA1, y = PCoA2, 
                     color = group_in_function),
                 level = 0.95, linetype = 2,linewidth = 0.8,
                 show.legend = FALSE) +
    geom_point(data = group_average, 
               aes(x = PCoA1, y = PCoA2, 
                   color = group_in_function),
               size = 2.2, show.legend = FALSE)+
    scale_x_continuous(limits = c(xmin,xmax))+
    scale_y_continuous(limits = c(ymin,ymax))+
    guides(color=guide_legend(override.aes = list(size=2.5)))+
    labs(x = xtitle,y = ytitle,color = '')+
    scale_color_manual(values = color)+
    annotate('text', label = label, 
             x = xmax*0.86, y = ymax*0.88, 
             size = 4.1, family = 'serif')+
    theme_bw()+
    theme(legend.key = element_blank(),
          legend.position = legend.position,
          legend.background = element_blank(),
          plot.background = element_blank(),
          # plot.title = element_text(family = "serif",size=15),
          axis.title.x= element_text(family = "serif",size=15),
          axis.title.y= element_text(family = "serif",size=15),
          axis.text.x= element_text(family = "serif",size=12),
          axis.text.y= element_text(family = "serif",size=12),
          legend.text = element_text(family = "serif",size=10),
          axis.ticks.length.x = unit(0.1,'cm'),
          axis.ticks.length.y = unit(0.1,'cm'),
          plot.margin = unit(c(0,0,0,0),"cm"))
  
  plot <- plot2+plot_spacer()+plot1+plot3+
    plot_layout(heights = c(1,4),widths = c(4,1),
                ncol = 2,nrow = 2,guides = 'keep')
  
  if(!is.null(path)){
    ggsave(plot = plot,
           filename = paste0(path,'/',filename,'_plot.pdf'),
           width = width,height = height)
  }
  
  return(plot)
}


