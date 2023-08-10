# Jialin He

\# function

\# 练习时长两年半，报错了必须是你的问题 o.0

\# 你怎么知道我抽到了卡芙卡！

## diversity-code

\# **metadata的row.names必须为sample ID**

\# **microdat的sample ID、OTU/Species/Genus...必须是行或列名**

\# color最好还是自定义

### α多样性

- **alpha_diversity**

```ruby
alpha_diversity(microdat = ,# 菌群特征表
                metadata = ,# 临床数据，row.names需要为sample ID
                group = ,# 分组变量
                plot_index = 'Shannon', # 默认'Shannon'，还有Simpson、Gini_Simpson、Inv_Simpson、Richness、Pielou、equitability、AVD
                p.adj = T,# T or F，是否校正P值，default = T
                p.signif = T,# T or F，default = T，T即返回***, **, *, ns，F返回P值
                title = NULL,# 图标题，可以为空
                color = NULL,# 可以为空
                path = NULL,# 保存pdf文件的路径，可以为空，即不输出到路径
                filename = 'diversity_plot',
                width = ,# plot宽度
                height = ) # plot高度
```


example:
```
res <- alpha_diversity(css_g,metadata,group = 'test',
                       color = color,plot_index = 'AVD',
                       p.adj = T,p.signif = T,
                       path = 'C://xx/')

# 返回sample ID、全部指数、分组的表格
res$data

# 返回ggplot图
res$plot
```
----

### β多样性

\# 非限制性，类似于无监督

- **PCOA_plot()**

属于经典多维排列 (Multidimensional scaling, MDS)分析

```ruby
PCOA_plot(microdat = ,# 菌群特征表
          metadata = ,# 临床数据，row.names需要为sample ID
          group = ,# 分组变量
          signif_method = , # adonis、anosim，可以为空
          distance = 'bray'
          title = NULL,# 图标题，可以为空
          color = NULL,# 可以为空
          legend.position = , # ggplot图例
          seed = 0,# 可以为空
          path = NULL,# 可以为空，保存pdf文件的路径，空即不输出到路径
          filename = 'PCOA',# 'PCOA'
          width = ,# plot宽度
          height = ) # plot高度
```

example:

```
res <- PCOA_plot(css_g,metadata,group = 'test',
                 color = color,path = 'C://xx/')

res$plot
```
![PCoA](https://github.com/JialinHe0o0/plot/blob/main/diversity_plot/PCoA_plot.png)

- **bray_dissimilarity()**

嘻嘻，你觉得相异距离只能用主成分展示吗

![bray-curtis](https://github.com/JialinHe0o0/plot/blob/main/diversity_plot/bray_curtis_plot.png)

- **PCoA_plus()**

  \# 你猜我放哪了

![](https://github.com/JialinHe0o0/plot/blob/main/diversity_plot/PCoA_plus_plot.png)

- **NMDS()**

非度量多维排列 (Non-metric multidimensional scaling, NMDS)
将距离矩阵转换为秩矩阵，更关注排序关系

```ruby
NMDS(microdat = ,
     metadata = ,
     group = ,
     distance = , # metaMDS()的参数，已有默认设置，可以为空
     k = ,try = ,trymax = , # metaMDS()的参数，已有默认设置，可以为空
     autotransform = , # metaMDS()的参数，已有默认设置，可以为空
     parallel = ,
     signif_method = , # adonis、anosim，可以为空
     title = ,
     color = ,
     legend.position = , # ggplot图例
     path = ,
     filename = ,
     seed = ,
     width = ,
     height = )
```

example

```ruby
res <- NMDS(css_g,metadata,group = 'test',
            signif_method = 'anosim',color = color,
            path = 'C://xx/')

res
```
![NMDS](https://github.com/JialinHe0o0/plot/blob/main/diversity_plot/NMDS_plot.png)


\# 限制性，类似于有监督

- **CPCOA_plot()**

```ruby
CPCOA_plot(microdat = ,# 菌群特征表
           metadata = ,# 临床数据，row.names需要为sample ID
           group = ,# 分组变量
           distance = 'bray'
           title = NULL,# 图标题，可以为空
           color = NULL,# 可以为空
           legend.position = , # ggplot图例
           seed = 0,# 可以为空
           path = NULL,# 可以为空，保存pdf文件的路径，空即不输出到路径
           filename = 'CPCOA',# 'CPCOA'
           width = ,# plot宽度
           height = ) # plot高度
```

example:

```ruby
res <- CPCOA_plot(css_g,metadata,group = 'test',
                  color = color,path = 'C://xx/')

res
```
![CPCoA](https://github.com/JialinHe0o0/plot/blob/main/diversity_plot/CPCoA_plot.png)

- **beta_plus()**

\# 基于betapart: https://doi.org/10.1111/j.2041-210X.2012.00224.x

\# 你可曾听过beta diversity的分解

\# 又是plus嘻嘻嘻

![turnover](https://github.com/JialinHe0o0/plot/blob/main/diversity_plot/beta1.png)
![nestedness](https://github.com/JialinHe0o0/plot/blob/main/diversity_plot/beta2.png)

