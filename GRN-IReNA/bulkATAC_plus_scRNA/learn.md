# bulkATAC和scRNA联合分析揭示基因调控网络
fimo是基于motif的PWM(position weight matrix)去扫描参考基因组数据，查看潜在的motif位置
将潜在的motif位点与bulkATAC的peaks取交集，则认为是真实稳定的motif，相对的target基因则是稳定的
然后使用cal_footprint_cut找到基于ATAC的potential_regulation，再将与基于scRNA用GENIE3推断的regulatory_relationships取交集得到最终的基因调控网络信息

[这个调控因子到底调不调控我的基因啊！——小果带你来学习FIMO数据库](https://mp.weixin.qq.com/s/_uYuLDY4_0E7tMQb4sGQNw)
[医学单细胞及表观多组学技术应用线上公开课(武汉大学/菲沙基因)整理总结及学习(第三部分-细胞互作/stripe/ATAC)](https://mp.weixin.qq.com/s/UbRLHbHiBTE5Lk6jEjx6XA)
[fimo扫描motif时的二三事](https://mp.weixin.qq.com/s/j3YSw-7Am4vwI321jwBQLQ)