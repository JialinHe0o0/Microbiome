
## detect_rate()

```ruby
detect_rate(dat = , # microdata丰度矩阵，可以是count或任何标准化方法处理后的data
            sample_in_row = T, # row.names == sample ID即为T
            proportion = 0.05) # 检出率，0.05即5%
```

## abun_filter()

注意，输入的丰度矩阵只能是count或者TSS标准化后的矩阵

```ruby
abun_filter(dat,
            sample_in_row = T,
            count = F, # 是否是count数据
            abundance = 0.001) # 平均丰度>0.001
```





