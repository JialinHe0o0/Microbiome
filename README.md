# Jialin He

\# function
\# 练习时长两年半

## diversity-code

\# **metadata的row.names必须为sample ID**

\# **microdat的sample ID、OTU/Species/Genus...必须是行或列名**

\# color最好还是自定义

### α多样性

- **alpha_diversity**

alpha_diversity(

microdat = ,# 菌群特征表

metadata = ,# 临床数据，row.names需要为sample ID

group = ,# 分组变量

sample_in_row = ,# T or F，T即特征表row.names为sample ID

plot_index = , # 默认'Shannon'，还有Simpson、Gini_Simpson、Inv_Simpson、Richness、Pielou、equitability、AVD

p.adj = ,# T or F，是否校正P值，default = T

p.signif = ,# T or F，default = T，T即返回***, **, *, ns，F返回P值

title = ,# 图标题，可以为空

color = ,# 可以为空

path = ,# NULL，保存pdf文件的路径，可以为空，即不输出到路径

filename = ,# 'shannon'

width = ,# plot宽度

height = ,# plot高度）

example:

res <- shannon_index(microdat = microdat, metadata = metadata, group = 'DM', sample_in_row = T, p.adj = T, title = 'Bacterial')

\# 返回sample ID、全部指数、分组的表格

res$data

\# 返回ggplot图

res$plot

----

### β多样性

\# 非限制性，类似于无监督

- **PCOA_plot()**

PCOA_plot(

microdat = ,# 菌群特征表

metadata = ,# 临床数据，row.names需要为sample ID

group = ,# 分组变量

sample_in_row = ,# T or F，T即特征表row.names为sample ID

title = ,# 图标题，可以为空

color = ,# 可以为空

seed = ,# 可以为空

path = ,# 可以为空，保存pdf文件的路径，空即不输出到路径

filename = ,# 'PCOA'

width = ,# plot宽度

height = ,# plot高度）


example:

res <- PCOA_plot(microdat, metadata, group = 'DM', sample_in_row = T)

res


\# 限制性，类似于有监督

- CPCOA_plot()

CPCOA_plot(

microdat = ,# 菌群特征表

metadata = ,# 临床数据，row.names需要为sample ID

group = ,# 分组变量

sample_in_row = ,# T or F，T即特征表row.names为sample ID

title = ,# 图标题，可以为空

color = ,# 可以为空

seed = ,# 可以为空

path = ,# 可以为空，保存pdf文件的路径，空即不输出到路径

filename = ,# 'CPCOA'

width = ,# plot宽度

height = ,# plot高度）

example:

res <- CPCOA_plot(microdat, metadata, group = 'DM', sample_in_row = T)

res


