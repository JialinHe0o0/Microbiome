
## detect_rate()

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

注意，输入的物种丰度矩阵只能是count或者TSS标准化后的矩阵

```ruby
abun_filter(dat, # 物种丰度矩阵
            sample_in_row = T, # T or F，T即特征表row.names为sample ID
            count = F, # 是否是count数据
            abundance = 0.001) # 平均丰度>0.001
```





