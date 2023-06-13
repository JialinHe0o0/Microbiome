# diversity-code

\# color最好还是自定义

\# metadata的row.names必须为sample ID

\# microdat的sample ID、OTU/Species/Genus...必须是行或列名

- α多样性

1. shannon_index()

> microdat = 菌群特征表

> metadata = 临床数据，row.names需要为sample ID

> group = 分组变量

> sample_in_row = T or F，特征表row.names为sample ID or colnames为sample ID

> p.adj = T or F，是否校正P值，default = T

> title = 图标题，可以为空

> color = ，可以为空

example:

res <- shannon_index(microdat = microdat, metadata = metadata, group = 'DM', sample_in_row = T, p.adj = T, title = 'Bacterial')

\# 返回sample ID与shannon index的表格

res$data

\# 返回ggplot图

res$plot




- β多样性

1. PCOA_plot()

> microdat = 菌群特征表

> metadata = 临床数据，row.names需要为sample ID

> group = 分组变量

> sample_in_row = T or F，特征表row.names为sample ID or colnames为sample ID

> title = 图标题，可以为空

> color = ，可以为空

> seed = ,可以为空



example:

res <- PCOA_plot(microdat, metadata, group = 'DM', sample_in_row = T)

res


