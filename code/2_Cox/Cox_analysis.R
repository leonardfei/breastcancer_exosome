aa=read.csv("riskscore_clinical.csv",header = T,row.names = 1)

library(survival)
library(plyr)
covariates <- colnames(aa)[c(4:11)]
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(X_OS, X_EVENT)~', x)))
univ_models <- lapply( univ_formulas, function(x){coxph(x, data = aa)})

#提取HR，95%置信区间和p值
univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         #获取p值
                         p.value<-signif(x$wald["pvalue"], digits=2)
                         #获取HR
                         HR <-signif(x$coef[2], digits=2);
                         #获取95%置信区间
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                         HR <- paste0(HR, " (", 
                                      HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(p.value,HR)
                         names(res)<-c("p.value","HR (95% CI for HR)")
                         return(res)
                       })
#转换成数据框，并转置
res <- t(as.data.frame(univ_results, check.names = FALSE))
as.data.frame(res)

#多因素cox
res.cox <- coxph(Surv(X_OS, X_EVENT) ~ riskscore + AJCC_Stage + Age + ER_Status + PR_Status + HER2_Final_Status+ Metastasis_Coded+ Node_Coded, data =  aa)
x <- summary(res.cox)
pvalue=signif(as.matrix(x$coefficients)[,5],2)
HR=signif(as.matrix(x$coefficients)[,2],2)
low=signif(x$conf.int[,3],2)
high=signif(x$conf.int[,4],2)
multi_res=data.frame(p.value=pvalue,
                     HR=paste(HR," (",low,"-",high,")",sep=""),
                     stringsAsFactors = F
)
multi_res
