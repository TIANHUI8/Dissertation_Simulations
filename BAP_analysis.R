setwd("/home/tianhui/GWAS/")
.libPaths(c('/home/tianhui/MyRlibs/',.libPaths()))
#########read data into R#####################################################################
###########################################################################################
Phenotype = read.csv(file = "BAP_phenotypes.csv", header = T, sep = ",", quote = "\"", dec = ".", fill = TRUE);
Markers = read.table(file = "BAP_impute_numeric.hmp.txt", header = T, stringsAsFactors = FALSE) ;
M = Markers[, -c(1:11)];
X = t(M);
XX = data.frame(rownames(X), X)
names(XX)[1] = 'Taxa'
Y1 = Phenotype[,-3];
mergedata1 = merge(XX, Y1, by = 'Taxa', all = TRUE, sort = FALSE)
newdata1 = mergedata1[1:345, ]
mydata1 = na.omit(newdata1)
X1 = mydata1[, 2:(ncol(mydata1)-1)]
Y1 = mydata1$height
#########single marker regression##############################################################
################################################################################################
p = ncol(X1);
pValues = numeric();
for(i in 1:p){
    fm = lm(Y1 ~ X1[, i])
    pValues[i] = summary(fm)$coef[2,4]
}
rank = order(pValues);
result1 = Markers[rank, 1]
####LASSO#########################################################################################
#####################################################################################################
library(glmnet)
set.seed(12356)

fit1 = cv.glmnet(as.matrix(X1), Y1, alpha = 1)
fit2 = glmnet(as.matrix(X1), Y1, alpha = 1, lambda = fit1$lambda.1se)
betahat = coef(fit2)[-1]
id = order(abs(betahat), decreasing = TRUE)
result2 = Markers[id, 1]
#####Elatic Net#####################################################################################
######################################################################################################
fit.en1 = cv.glmnet(as.matrix(X1), Y1, alpha = 0.5)
fit.en2 = glmnet(as.matrix(X1), Y1, alpha = 0.5, lambda = fit.en1$lambda.1se)
betahat2 = coef(fit.en2)[-1]
id2 = order(abs(betahat2), decreasing = TRUE)
result3 = Markers[id2, 1]
####Ridge Regression###############################################################################
#######################################################################################################
fit.RR1 = cv.glmnet(as.matrix(X1), Y1, alpha = 0)
fit.RR2 = glmnet(as.matrix(X1), Y1, alpha = 0,lambda = fit.RR1$lambda.1se)
betahat3 = coef(fit.RR2)[-1]
id3 = order(abs(betahat3), decreasing = TRUE)
result4 = Markers[id3, 1]
########Bayesian Alphabet#############################################################################
###########################################################################################################
nIter = 12000; burnIn = 2000
library(BGLR)
set.seed(12345)
X1 = scale(X1, center = TRUE, scale = TRUE)
####Bayesian Ridge regression#############
ETA = list(MRK = list(X = X1, model = 'BRR'))
fmBRR = BGLR(y = Y1, ETA = ETA, nIter = nIter, burnIn = burnIn, saveAt = 'BRR_')

a = fmBRR$ETA[[1]]$b
id5 = order(abs(a), decreasing = TRUE)
result5 = Markers[id5, 1]
###BayesA(Scaled-t prior)###
ETA$MRK$model = 'BayesA'
fmBA = BGLR(y = Y1, ETA = ETA, nIter = nIter, burnIn = burnIn, saveAt = 'BA_')

b = fmBA$ETA[[1]]$b
id6 = order(abs(b), decreasing = TRUE)
result6 = Markers[id6, 1]

result = cbind(result1, result2, result3, result4, result5, result6)
save(result, file = "BAP_result.RData")

