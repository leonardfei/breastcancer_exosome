suppressMessages(library(GSVA))
suppressMessages(library(GSVAdata))
suppressMessages(library(GSEABase))
suppressMessages(library(limma))

#读入gmt文件，这个可以从MSigDB上下载，这边选的上gene symbol根据自己的data来选择
gmt_file="~/Downloads/c5.go.v7.4.symbols.gmt"
geneset <- getGmt(gmt_file)  
#自行读入exp文件，我这边为行为gene symbol，列为样本
es <- gsva(as.matrix(exp), geneset,
           min.sz=10, max.sz=500, verbose=TRUE,kcdf="Poisson")#如果是芯片数据，log2即可用；如果是转录组数据，则count，RPKM，FPKM，TPM皆可用，使用count时须指明参数kcdf=“Poisson”
#得到gsva计算的数值后再用limma包做差异分析得到差异的pathway
design <- model.matrix(~ factor(c(rep("cont",4),rep("treatment",5))))
colnames(design) <- c("ALL", "contvstreatment")
row.names(design)<-colnames(exp)
fit <- lmFit(es, design)
fit <- eBayes(fit)
#这边是总的
allGeneSets <- topTable(fit, coef="contvstreatment", number=Inf)
#这边是差异的
adjPvalueCutoff <- 0.001
DEgeneSets <- topTable(fit, coef="contvstreatment", number=Inf,
                       p.value=adjPvalueCutoff, adjust="BH")
res <- decideTests(fit, p.value=adjPvalueCutoff)

####################GSEA################
library(edgeR)
library(clusterProfiler)

exprSet=exp
group_list=c(rep("normal",number of low_group),rep("cancer",number of high_group))
d <- DGEList(counts=exprSet,group=factor(group_list))#构建模型，只需要两个因素
keep <- rowSums(cpm(d)>1) >= 2#做筛选，此方法需要这样做
table(keep)
d <- d[keep, , keep.lib.sizes=FALSE]
d$samples$lib.size <- colSums(d$counts)
d <- calcNormFactors(d)#TMM归一化
d$samples
dge=d
design <- model.matrix(~0+factor(group_list))
rownames(design)<-colnames(dge)
colnames(design)<-levels(factor(group_list))
dge=d 
dge <- estimateGLMCommonDisp(dge,design)
dge <- estimateGLMTrendedDisp(dge, design)
dge <- estimateGLMTagwiseDisp(dge, design)
fit <- glmFit(dge, design)  
# https://www.biostars.org/p/110861/
can.vs.nor <- makeContrasts(cancer-normal, levels=design)
lrt <- glmLRT(fit,contrast=can.vs.nor) #因为是两组所以这么做，按照说明书
nrDEG=topTags(lrt, n=nrow(dge))
nrDEG=as.data.frame(nrDEG)
head(nrDEG)

gene=nrDEG
geneList=gene$logFC
names(geneList)=rownames(gene)
geneList = sort(geneList, decreasing = TRUE)
kegmt<-read.gmt("h.all.v2023.1.Hs.symbols.gmt")
KEGG<-GSEA(geneList,TERM2GENE = kegmt,eps = 0,seed = 123,pvalueCutoff = 1)

library(enrichplot)
gseaplot2(KEGG,1,color="red",pvalue_table = T)
out=KEGG@result
write.csv(out,"result.csv")




