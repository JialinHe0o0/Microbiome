# Jialin He

## Overview

- **特征与群落变化的关联**
- **特征对群落变化的解释度**
- **特征重要性分解**

## mantel_function()

```ruby
mantel_function(microdat, # 菌群特征表，也就是待分析的群落
                metadata, # 临床数据，row.names需要为sample ID
                sample_in_row=T, # T or F，T即菌群特征表row.names为sample ID
                var, # 特征，例如c('age','BMI','ABSI')
                method = 'pearson', # 'pearman' or 'spearman'
                community_name='Fungi', # 群落名称
                permutations = 999,
                parallel = 4, # 线程
                seed = 0,
                rbreak = c(-Inf,0.1,0.2,Inf), # 相关系数切点，方便后续作为图形的legend
                rlabel = c("<0.1","0.1-0.2",">=0.2"),
                pbreak = c(-Inf,0.001,0.01,0.05,Inf),  # P值切点，方便后续作为图形的legend
                plabel = c("<0.001","<0.01","<0.05",">0.05"))
```

rbreak, rlabel, pbreak, plabel并不重要，只是后续作图时，可以用切割后的r, P值来规定曲线颜色、类型等，方便展示

完全可以出结果后修改

example

```ruby
var <- c('age','BMI','waist','TG','TC','HDL','LDL')

res <- mantel_function(microdat = asv_css,metadata = table1,var = var,
                       method = 'pearson',sample_in_row = T,
                       community_name = 'Fungi',
                       rbreak = c(-Inf,0.05,0.1,Inf),
                       rlabel = c('<0.05','0.05-0.1','>0.1'))

res2 <- mantel_function(microdat = bac_species,
                        metadata = table1,var = var,
                        method = 'pearson',sample_in_row = T,
                        community_name = 'Bacteria',
                        rbreak = c(-Inf,0.05,0.1,Inf),
                        rlabel = c('<0.05','0.05-0.1','>0.1'))

res_all <- rbind(res,res2)

# 绘图
library(linkET)
# 参考: https://github.com/Hy4m/linkET
# 嘿嘿，不能不劳而获

```

最终绘图

（参数自有道理）

![Mantel_test](https://github.com/JialinHe0o0/Microbiome/blob/main/plot/Mantel_test.png)

怎么解释呢，嘿嘿，自己想





