# Jialin He

\# function

\# 练习时长两年半，报错了必须是你的问题 o.0

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
                sample_in_row = T,# T or F，T即特征表row.names为sample ID
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
res <- alpha_diversity(microdat = bacteria,
                       metadata = metadata,
                       group = 'Group',
                       sample_in_row = T,p.adj = T,
                       title = 'Example',p.signif = T,
                       color = color[1:5])

# 返回sample ID、全部指数、分组的表格
res$data

# 返回ggplot图
res$plot
```

----

### β多样性

\# 非限制性，类似于无监督

- **PCOA_plot()**

```ruby
PCOA_plot(microdat = ,# 菌群特征表
          metadata = ,# 临床数据，row.names需要为sample ID
          group = ,# 分组变量
          sample_in_row = T,# T or F，T即特征表row.names为sample ID
          title = NULL,# 图标题，可以为空
          color = NULL,# 可以为空
          seed = 0,# 可以为空
          path = NULL,# 可以为空，保存pdf文件的路径，空即不输出到路径
          filename = 'PCOA',# 'PCOA'
          width = ,# plot宽度
          height = ) # plot高度
```

example:

```
res <- PCOA_plot(microdat, metadata, group = 'DM', sample_in_row = T)

res
```

- NMDS()

```ruby
NMDS(microdat = ,
     metadata = ,
     group = ,
     sample_in_row = ,
     title = ,
     color = ,
     path = ,
     distance = , # metaMDS()的参数，已有默认设置，可以为空
     k = ,try = ,trymax = , # metaMDS()的参数，已有默认设置，可以为空
     autotransform = , # metaMDS()的参数，已有默认设置，可以为空
     parallel = ,
     filename = ,
     seed = ,
     width = ,
     height = )
```

example

```ruby
res <- NMDS(microdat, metadata, group = 'DM', sample_in_row = T)
```


\# 限制性，类似于有监督

- CPCOA_plot()

```ruby
CPCOA_plot(microdat = ,# 菌群特征表
           metadata = ,# 临床数据，row.names需要为sample ID
           group = ,# 分组变量
           sample_in_row = T,# T or F，T即特征表row.names为sample ID
           title = NULL,# 图标题，可以为空
           color = NULL,# 可以为空
           seed = 0,# 可以为空
           path = NULL,# 可以为空，保存pdf文件的路径，空即不输出到路径
           filename = 'CPCOA',# 'CPCOA'
           width = ,# plot宽度
           height = ) # plot高度
```

example:

```
res <- CPCOA_plot(microdat, metadata, group = 'DM', sample_in_row = T)

res
```



