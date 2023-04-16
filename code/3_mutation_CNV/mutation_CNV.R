BiocManager::install("maftools")

library(maftools)
#设置路径获取R包内置的 TCGA LAML MAF 文件
#laml.maf = system.file('extdata', 'tcga_laml.maf.gz', package = 'maftools')
#clinical information containing survival information and histology. This is optional
#laml.clin = system.file('extdata', 'tcga_laml_annot.tsv', package = 'maftools')
#laml = read.maf(maf = laml.maf, clinicalData = laml.clin)
#读取maf文件进行个性化处理maf=read.table("TCGA.BRCA.varscan.6c93f518-1956-4435-9806-37185266d248.DR-10.0.somatic.maf",header = T,sep = "\t",fileEncoding = "GBK",fill = T,quote = "")
#写出来maf文件write.table(maf_low,"maf_low.maf",sep = "\t",fileEncoding = "GBK",row.names = F,col.names = T,quote = F)

########################low risk################3
laml.maf="/Users/mac/Desktop/数据库分析/深度分析/突变/gdc_download_20210705_014158.760126/maf_low.maf"
laml = read.maf(maf = laml.maf)
x = tmb(maf = laml)#突变负荷
plotmafSummary(maf = laml, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)

#oncoplot for top ten mutated genes.
oncoplot(maf = laml, top = 10,borderCol=NULL)
#oncostrip(maf = laml, genes = c('DNMT3A','NPM1', 'RUNX1'))
#oncoplot(maf = laml, genes = c("S100A1","S100A2","S100A3","S100A4","S100A5","S100A6","S100A7","S100A8","S100A9","S100A10","S100A11","S100A12","S100A13","S100A14","S100A16","S100A7A","S100A7L2","S100B","S100G","S100P","S100Z"))

laml.titv = titv(maf = laml, plot = FALSE, useSyn = TRUE)
#plot titv summary
plotTiTv(res = laml.titv)
laml_L=laml
########################high risk################3
laml.maf="/Users/mac/Desktop/数据库分析/深度分析/突变/gdc_download_20210705_014158.760126/maf_high.maf"
laml = read.maf(maf = laml.maf)
x = tmb(maf = laml)#突变负荷
plotmafSummary(maf = laml, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)

#oncoplot for top ten mutated genes.
oncoplot(maf = laml, top = 10,borderCol=NULL)
#oncostrip(maf = laml, genes = c('DNMT3A','NPM1', 'RUNX1'))
#oncoplot(maf = laml, genes = c("S100A1","S100A2","S100A3","S100A4","S100A5","S100A6","S100A7","S100A8","S100A9","S100A10","S100A11","S100A12","S100A13","S100A14","S100A16","S100A7A","S100A7L2","S100B","S100G","S100P","S100Z"))

laml.titv = titv(maf = laml, plot = FALSE, useSyn = TRUE)
#plot titv summary
plotTiTv(res = laml.titv)
laml_H=laml

#lollipop plot for DNMT3A, which is one of the most frequent mutated gene in Leukemia.
#lollipopPlot(maf = laml, gene = 'DNMT3A', AACol = 'HGVSp_Short', showMutationRate = TRUE)
lollipopPlot2(m1 = laml_H, m2 = laml_L, gene = "TP53",AACol1 = "HGVSp_Short", AACol2 = "HGVSp_Short", m1_name = "exoSIG High", m2_name = "exoSIG Low",domainLabelSize = 1,labPosSize=10,legendTxtSize=1)

#gistic结果文件的CNV分析

laml.gistic = readGistic(gisticAllLesionsFile ='/Users/mac/Desktop/数据库分析/深度分析/CNV/gistic2_outdir_CNV_high/all_lesions.conf_90.txt' , gisticAmpGenesFile ='/Users/mac/Desktop/数据库分析/深度分析/CNV/gistic2_outdir_CNV_high/amp_genes.conf_90.txt' , gisticDelGenesFile ='/Users/mac/Desktop/数据库分析/深度分析/CNV/gistic2_outdir_CNV_high/del_genes.conf_90.txt' , gisticScoresFile ='/Users/mac/Desktop/数据库分析/深度分析/CNV/gistic2_outdir_CNV_high/scores.gistic' , isTCGA = TRUE)
#> -Processing Gistic files..
#> --Processing amp_genes.conf_99.txt
#> --Processing del_genes.conf_99.txt
#> --Processing scores.gistic
#> --Summarizing by samples
#GISTIC object
laml.gistic

gisticChromPlot(gistic = laml.gistic, markBands = "all",ref.build = "hg38")
gisticBubblePlot(gistic = laml.gistic)



