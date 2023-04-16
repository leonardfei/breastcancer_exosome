library(glmnet)
dev=read.csv("origin.csv",header = T,row.names = 1)
table(dev$type)

dev$type=c(rep(1,1049),rep(0,76))
x=dev[,1:16]
y=dev[,17]
x=as.matrix(x)
fit<-glmnet(x,y,alpha=1,family='binomial')
plot(fit, xvar = "lambda", label = TRUE)
print(fit)

cv.fit <- cv.glmnet(x,y,alpha=1,nlambda = 1000)
plot(cv.fit)
abline(v=log(c(cv.fit$lambda.min,cv.fit$lambda.1se)),lty=2)

#如果取最小值时
cv.fit$lambda.min
Coefficients <- coef(fit, s = cv.fit$lambda.min)
Active.Index <- which(Coefficients != 0)
Active.Coefficients <- Coefficients[Active.Index]
Active.Index
Active.Coefficients
row.names(Coefficients)[Active.Index]


#如果取1倍标准误
cv.fit$lambda.1se
Coefficients <- coef(fit, s = cv.fit$lambda.1se)
Active.Index <- which(Coefficients != 0)
Active.Coefficients <- Coefficients[Active.Index]
Active.Index
Active.Coefficients
row.names(Coefficients)[Active.Index]