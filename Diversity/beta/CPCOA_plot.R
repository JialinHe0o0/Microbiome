# Jialin He 202308

CPCOA_plot <- function(microdat,metadata,group,
                       distance = 'bray',
                       parallel = 4,
                       title = NULL,
                       color = NULL,
                       legend.position = 'top',
                       path = NULL,
                       filename = 'CPCoA',
                       seed = 0,
                       width = 6.6,
                       height = 5.1){
  
  if(!require(pacman))install.packages(pacman)
  pacman::p_load(tidyverse,vegan,car,ggthemes,scico)

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
  
  # distance matrix
  dis_mat <- vegdist(microdat, method = distance) %>% as.matrix()
  
  len <- length(unique(metadata$group_in_function))
  
  if(is.null(color)){
    color <- scico::scico(n = len,palette = 'berlin')
  }
  
  if (len < 3){
    stop(print('ERROR: You need at least 3 levels in the group'))
  }else{
    # 函数提取CCA中主要结果
    variability_table = function(cca){
      chi = c(cca$tot.chi, cca$CCA$tot.chi, cca$CA$tot.chi)
      variability_table = cbind(chi, chi/chi[1])
      colnames(variability_table) = c("inertia", "proportion")
      rownames(variability_table) = c("total", "constrained", "unconstrained")
      return(variability_table)
    }
    
    # Constrained analysis OTU table by genotype
    capscale.gen = capscale(as.dist(dis_mat) ~ group_in_function, 
                            data=metadata, add=F, sqrt.dist=T)
    
    # ANOVA-like permutation analysis, calculate P value
    set.seed(seed)
    perm_anova.gen = anova.cca(capscale.gen, permutations = 1000, 
                               parallel = parallel)
    p.val = round(perm_anova.gen[1,4],3)
    
    # generate variability tables and calculate confidence intervals for the variance
    var_tbl.gen = variability_table(capscale.gen)
    eig = capscale.gen$CCA$eig
    
    # 解释度
    variance = round(var_tbl.gen["constrained", "proportion"]*100,3)

    
    # extract the weighted average (sample) scores
    points = as.data.frame(capscale.gen$CCA$wa)
    points = merge(metadata,points,'row.names')
    points <- points[,c('Row.names','group_in_function','CAP1','CAP2')]
    
    group_average <- aggregate(cbind(CAP1, CAP2)~group_in_function, 
                               data = points, FUN = mean)
    
    xmin <- NULL
    xmax <- NULL
    ymin <- NULL
    ymax <- NULL
    
    for(i in unique(metadata$group_in_function)){
      temp <- points[points$group_in_function == i,]
      eli <- dataEllipse(x = temp$CAP1,y = temp$CAP2,
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
    
    variance <- paste0(variance,'%')
    label <- bquote(atop(R^2 == .(variance),  italic(P) ==.(p.val)))
    
    # plot CPCoA 1 and 2
    plot <- ggplot(points, aes(x = CAP1, y = CAP2, 
                           color = group_in_function)) + 
      geom_point(size = 0.8) +
      labs(x=paste("CPCoA1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
           y=paste("CPCoA2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) +
      stat_ellipse(level = 0.95, linetype = 2,linewidth = 0.8,
                   show.legend = FALSE) +
      geom_point(data = group_average, 
                 aes(x = CAP1, y = CAP2, 
                     color = group_in_function),
                 size = 2.2, show.legend = FALSE)+
      scale_x_continuous(limits = c(xmin-(xmax-xmin)*0.025,
                                    xmax+(xmax-xmin)*0.025))+
      scale_y_continuous(limits = c(ymin-(ymax-ymin)*0.025,
                                    ymax+(ymax-ymin)*0.025))+
      guides(color=guide_legend(override.aes = list(size=2.5)))+
      labs(title = title,color = '')+
      scale_color_manual(values = color)+
      # annotate('text', label = 'PERMANOVA', 
      #          x = xmax*0.9, y = xmax*0.99, size = 3.6)+
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
  }
  
  if(!is.null(path)){
    ggsave(plot = plot,
           filename = paste0(path,'/',filename,'_plot.pdf'),
           width = width,height = height)
  }
  
  return(plot)
}
