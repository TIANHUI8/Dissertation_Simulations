setwd("D:/Clemson/RESEARCH//RESEARCH/mypaper/simulation/all243results")
library(glmnet)
library(BGLR)
sims = 1000   #number of replication
p = 1000  # number of markers
n = 100   #number of individual
nIter = 12000; burnIn = 2000
q = runif(p, min = 0, max = 1)
cutoff1 = seq(5*10^(-5), 0.05, by = 10^(-5))
#cutoff2 = seq(0, 1, by = 0.01)
cutoff3 = seq(0.01, 1, by = 0.01)
len1 = length(cutoff1)
#len2 = length(cutoff2)
len3 = length(cutoff3)
p.cutoff = matrix(0, ncol = len1, nrow = sims)
#alpha.cutoff = matrix(0, ncol = len2, nrow = sims)
beta3.cutoff = matrix(0, ncol = len3, nrow = sims)
for(s in 1:sims){
        X = t(matrix(rbinom(n*p, 2, rep(q, n)), ncol = n))
        Y =  3 + rnorm(n, 0, 1) 
        ###single marker regression##
        pValues=c(); beta.hat1 = c(); 
        for(a in 1:p){
             fm = lm(Y ~ X[,a])
	       if(length(unique(X[,a])) == 1){
                          pValues = c(pValues, 1)
                          beta.hat1 = c(beta.hat1, 0)
	       }
	       if(length(unique(X[,a])) > 1){
               pValues = c(pValues, summary(fm)$coef[2,4])
               beta.hat1 = c(beta.hat1, summary(fm)$coef[2,1])
	       }							
         }
         for(i in 1:len1){
               p.cutoff[s,i] = quantile(pValues, cutoff1[i])
               #s1[s,i] = sum(pValues < p.cutoff)   
         }                   
         ##########BayesA########################################            
         ETA = list(MRK = list(X = X, model = 'BayesA'))
         fmBA = BGLR(y = Y, ETA = ETA, nIter = nIter, burnIn = burnIn, saveAt = 'BA_')
         beta.hat3 = fmBA$ETA[[1]]$b
         for(i in 1:len3){
         beta3.cutoff[s,i] = quantile(abs(beta.hat3), 1 - cutoff3[i])
         }
}
write.csv(p.cutoff, file = "aoc_cutoff_pvalue1.csv")
write.csv(beta3.cutoff, file = "aoc_cutoff_beta.csv")