library(dplyr)
library(tools)
library(optparse)

## 1. 定义选项
option_list <- list(
  make_option(c("-d", "--dir_path"),
              type = "character",
              default = "/data/work/Single-Cell-Pipeline/GRN-gene2role/cotton/input",
              help = "directory containing *_el_signed.txt files [default= %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))
dir_path <- opt$dir_path

files <- list.files(dir_path, full.names = TRUE)

## 2. 循环读入，收集各文件基因集合
gene_sets <- list()
group <- c()
for (f in files) {
  dat <- read.table(f, header = FALSE, sep = "\t")
  gene_set <- unique(c(dat$V1, dat$V2))
  sample_name <- file_path_sans_ext(basename(f))  # 去掉扩展名
  gene_sets[[sample_name]] <- gene_set
  group <- c(group, sample_name)
}
# head(dat)
# #           V1            V2 V3
# # 1 Ga13g00904    Ga07g01276  1
# # 2 Ga04g00539    Ga06g00349  1
# # 3 Ga13g00777    Ga14g01413  1
# # 4 Ga04g00539    Ga10g02684  1
# # 5 Ga04g00539 lnc-Ga09g0097  1
# # 6 Ga04g00539    Ga04g01599  1

# ## 3. 求所有样本基因集合的交集
# common_genes <- Reduce(intersect, gene_sets) |> sort()

# 所有样本基因集合的并集（去重并按字母排序）
common_genes <- Reduce(union, gene_sets) |> sort()

## 4. 生成模板
# 已知变量
# common_genes <- c("Ga04g00539", "Ga07g00382", "Ga12g00591")
# group <- c("N2", "N1", "P0", "P1")

# 构造数值矩阵
n <- length(common_genes)
values <- matrix(1:(n * length(group)), nrow = n, byrow = FALSE)

# 构造data.frame
df <- data.frame(
  gene = common_genes,
  values
)

# 设置列名
# group <- c("day-1", "day-2", "day0", "day1")
colnames(df) <- c("gene", group)

# 查看结果
print(head(df))

## 1. 先把 df 变成“基因 → 数字”的查表向量
gene_map <- setNames(df[ , -1, drop = FALSE], group)   # 每个分组一列
gene_map <- lapply(gene_map, function(x) setNames(x, df$gene))

## 2. 构造一个函数：根据分组名和基因名返回对应数字
get_num <- function(vec, grp)
  unlist(gene_map[[grp]][ vec ])


## 3. 过滤 + 替换
result_all <- NULL   # 初始化为空

for (i in seq_along(files)) {
  dat <- read.table(files[i], header = FALSE, sep = "\t")
  tmp <- dat %>%
    filter(V1 %in% df$gene & V2 %in% df$gene) %>%
    mutate(V1_num = get_num(V1, group[i]),
           V2_num = get_num(V2, group[i]),
           V3_num = V3) %>%
    select(V1_num, V2_num, V3_num)
  result_all <- rbind(result_all, tmp)   # 逐次拼接
}


head(result_all)
write.table(result_all,file = "sample_merged.edgelist",sep = "\t",row.names = FALSE,col.names = FALSE,quote = FALSE)
                   
rownames(df) <- df$gene
df <- df[, -1]
## 5. 保存
# 如果 splitMatrix 不存在就创建
dir.create("splitMatrix", showWarnings = FALSE, recursive = TRUE)
write.table(df,file = "./splitMatrix/index_tracker.tsv",sep = "\t",row.names = TRUE,col.names = TRUE,quote = FALSE)