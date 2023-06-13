# diversity-code
- shannon_index()：α多样性
> microdat = 菌群特征表

> metadata = 临床数据，row.names需要为sample ID

> group = 分组变量

> sample_in_row = T or F，特征表row.names为sample ID or colnames为sample ID

> p.adj = T or F，是否校正P值

> title = 表格标题，可以为空

example:
res <- shannon_index(microdat = microdat, metadata = metadata, group = 'DM', sample_in_row = T, p.adj = T, title = 'Bacterial')

