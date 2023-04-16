#下面就是进行分析了，在进行分析的时候我们需要准备三个文件：

#LM22.txt（可以从CIBERSORT网站下载，这个就是22种免疫细胞的参考marker基因表达）
#CIBERSORT.R（CIBERSORT源代码，从github下载:https://rdrr.io/github/singha53/amritr/src/R/supportFunc_cibersort.R）
#expression.txt（基因表达谱文件）
#第一列是基因名，第一行是样品名，不能有重复基因名，第一列列名不能空白。矩阵中不能存在空白或NA值，不要对表达量取Log2.
# 1.不可以有负值和缺失值 2.不要取log 3.如果是芯片数据，昂飞芯片使用RMA标准化，Illumina 的Beadchip 和Agilent的单色芯片，用limma处理。 4.如果是RNA-Seq表达量，使用FPKM和TPM都很合适。
#如果表达矩阵中基因不能完全覆盖LM22.txt中的基因，Cibersort同样可以正常运行，但不能少于LM22.txt中所需基因的一半。
#表达矩阵保存为制表符分割的txt文本（“DATA.txt”）

#然后，运行CIBERSORT非常简单，三个文件放到一个文件夹中，在RStudio中设置工作路径，运行如下代码即可：
setwd("/Users/mac/Downloads/cibersort免疫浸润分析")
source("CIBERSORT.R")

# Define LM22 file
LM22.file <- "LM22.txt"
exp.file <- "exp.txt"

TME.results = CIBERSORT(LM22.file, exp.file, perm = 1000, QN = TRUE)#perm推荐1000，但时间变长，测序QN推荐F，芯片推荐T

# output CIBERSORT results
write.table(TME.results, "TME.results.output.txt", 
            sep = "\t", row.names = T, col.names = T, quote = F)

plot.info=TME.results

# boxplot
library(ggpubr)
library(ggsci)
plot.info=read.csv("plotinfo.csv",header = T)#该文件需要对TME结果进行处理成一列细胞类型，另一列表达值的格式
p=ggboxplot(
  plot.info,
  x = "CellType",
  y = "Composition",
  color = "CellType",
  xlab = "",
  ylab = "Cell composition",
  main = "TME Cell composition",
  add  = "jitter"
)  +
  theme(axis.text.x = element_text(
    angle = 90,
    hjust = 1,
    vjust = 1
  ))
p+scale_color_jco()


# Barplot of cell componment of each sample using hclust
TME.data=read.table("TME.results.output.txt",header = T,row.names = 1,sep = "\t")
sample.index <-
  hclust(dist(TME.data), method = "ward.D")$order
sample.order <- rownames(TME.data)[sample.index]

ggbarplot(
  plot.info,
  x = "PATIENT_ID",
  y = "Composition",
  size = 0,
  fill = "CellType",
  color = "CellType",
  order = sample.order
)  +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 1,
      size = 1
    ),
    legend.position = "bottom"
  )

#################ssGSEA####################
exp=read.table("exp.txt",header = T,row.names = 1,sep = "\t")
exp=2^exp-1
exp=as.matrix(exp)
#如果是芯片数据，log2即可用；如果是转录组数据，则count，RPKM，FPKM，TPM皆可用，使用count时须指明参数kcdf=“Poisson”。
library(GSVA)
re <- gsva(exp, geneset, method="ssgsea",
           mx.diff=FALSE, verbose=FALSE,kcdf="Poisson")
library(pheatmap)
re2 <- as.data.frame(re)

an = data.frame(group = Group,
                row.names = colnames(exp))
pheatmap(re2,scale = "row",
         show_colnames = F,
         #annotation_col = an,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100))

library(Hmisc)
identical(colnames(re),colnames(exp))
gs = c("RACGAP1","ADRB1","FLT3LG")
nc = t(rbind(re,exp[gs,]))
nc[1:4,1:4]

m = rcorr(nc)$r[1:nrow(re),(ncol(nc)-length(gs)+1):ncol(nc)]
p = rcorr(nc)$P[1:nrow(re),(ncol(nc)-length(gs)+1):ncol(nc)]

library(dplyr)
tmp = matrix(case_when(p<0.01~"**",
                       p<0.05~"*",
                       T~""),nrow = nrow(p))
tmp2=as.data.frame(t(tmp))
n=as.data.frame(t(m))

source("modified_pheatmap.R")
pheatmap(n,
         display_numbers =tmp2,
         angle_col =45,
         color = colorRampPalette(c("#92b7d1", "white", "#d71e22"))(100),
         border_color = "white",
         treeheight_col = 0,
         treeheight_row = 0)

