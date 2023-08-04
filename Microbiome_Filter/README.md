# Jialin He

# 高丰度/高检出率物种的筛选

## detect_rate()

基于物种的检出率得到core feature

```ruby
detect_rate(dat = , # 物种丰度矩阵，可以是count或任何标准化方法处理后的data
            sample_in_row = T, # T or F，T即特征表row.names为sample ID
            proportion = 0.05) # 检出率，0.05即5%
```

example

```ruby
res <- detect_rate(dat = css_g, sample_in_row = T, proportion = 0.05)

# 返回每个物种的检出率
res$rate

# 返回剔除了检出率低于5%的物种的矩阵
res$dat

```

## abun_filter()

基于物种的平均相对丰度得到core feature

注意，输入的物种丰度矩阵只能是count或者TSS标准化后的矩阵

```ruby
abun_filter(dat, # 物种丰度矩阵
            sample_in_row = T, # T or F，T即特征表row.names为sample ID
            count = F, # 是否是count数据
            abundance = 0.001) # 平均丰度>=0.001
```

example

```ruby
res <- abun_filter(tss_g, sample_in_row = T, count = F, abundance = 0.001)

# res即为剔除了平均相对丰度低于0.001的物种的矩阵

```

# 差异分析

一般而言，筛出高丰度/高检出率物种后，做组间差异分析

## DA_test()

最简单的差异分析方法

2组： t.test或wilcox.test

3组或以上: ANOVA或Kruskal.Wallis

```ruby
DA_test(microdata = , # 物种丰度矩阵，必须是标准化后的data
        metadata = , # clinical data
        group = , # 分组变量名称
        parametric_test=F，# 是否非参数检验
        cutoff_padj=0.05) # 校正后P值的切点
```

example

```ruby
res <- DA_test(microdat = dat_filter,
               metadata = metadata,
               group = 'Group',
               parametric_test = F，
               cutoff_padj=0.05)
res
```
|taxo|pvalue|p_adj|
|---|---|---|
|a菌|0.003|0.006|
|b菌|0.009|0.009|







