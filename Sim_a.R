library(glmnet)
library(BGLR)
sims = 4   #number of replication
p = 1000  # number of markers
n = c(10, 45, 100)  #number of individual
lenn = length(n)
h = c(200, 100, 0)  #number of misimputated marker
lenh = length(h)
nIter = 12000; burnIn = 2000
beta.true = c(0.25, 0.5, 0.75, 1, 2, 5, rep(0, p-6))
p.cutoff = c(1:p)*0.05/p
#beta.cutoff = 0.01563455
beta.cutoff = 0.03086782 
r = c(0.8, 0.4)
lenr = length(r)
fp1 = fp2 = fp3 = matrix(0, sims, 243)
tp1 = tp2 = tp3 = matrix(0, sims, 243)
fdr1 = fdr2 = fdr3 = matrix(0, sims, 243)
fpr1 = fpr2 = fpr3 = matrix(0, sims, 243)
tpr1 = tpr2 = tpr3 = matrix(0, sims, 243)
mse1 = mse2 = mse3 = matrix(0, sims, 243)
nmv1 = nmv2 = nmv3 = matrix(0, sims, 243)
s2 = s3 = matrix(0, sims, 243)
q = runif(p, min = 0, max = 1)
for(s in 1:sims){
count = 1
    for(i in 1:lenn){
        for(j in 1:lenh){
            for(l in 1:3){ # loop over cor_x
                 if(l == 3){
                    X = t(matrix(rbinom(n[i]*p, 2, rep(q, n[i])), ncol = n[i]))
                 }
                 else{
                    Z1 = rnorm(n[i], 0, 1)
                    Z2 = r[l]*Z1 + rnorm(n[i], 0, sd = sqrt(1 - r[l]^2))
                    Z.1 = c()
                    for(k in 1:n[i]){
	                  if(Z1[k] < -0.67){Z.1[k] = 0}  
	                  else if(Z1[k] >= -0.67 & Z1[k] <= 0.67){Z.1[k] = 1} 
                        else{Z.1[k] <- 2} 
                    }
                    Z.2 = c()
                    for(k in 1:n[i]){
                    if(Z2[k] < -0.67){Z.2[k] <- 0}  
                    else if(Z2[k] >= -0.67 & Z2[k] <= 0.67){Z.2[k] = 1} 
                    else{Z.2[k] <- 2} 
                    }
                    Z3 = rnorm(n[i], 0, 1)
                    Z4 = r[l]*Z3 + rnorm(n[i], 0, sd = sqrt(1-r[l]^2))
                    Z.3 = c()
                    for(k in 1:n[i]){
                        if(Z3[k] < -0.67){Z.3[k] <- 0}  
                        else if(Z3[k] >= -0.67 & Z3[k] <= 0.67){Z.3[k] <- 1} 
                        else{Z.3[k] <- 2} 
                    }
                    Z.4 = c()
                    for(k in 1:n[i]){
	                   if(Z4[k] < -0.67){Z.4[k] <- 0}  
	                   else if(Z4[k] >= -0.67 & Z4[k] <= 0.67){Z.4[k] <- 1} 
                         else{Z.4[k] <- 2} 
                    }
                    Z5 = rnorm(n[i], 0, 1)
                    Z6 = r[l]*Z5 + rnorm(n[i], 0, sd = sqrt(1-r[l]^2))
                    Z.5 = c()
                    for(k in 1:n[i]){
	              if(Z5[k] < -0.67){Z.5[k] <- 0}  
                    else if(Z5[k] >= -0.67 & Z5[k] <= 0.67){Z.5[k] <- 1} 
                    else{Z.5[k] <- 2} 
                    }
                    Z.6 = c()
                    for(k in 1:n[i]){
	              if(Z6[k] < -0.67){Z.6[k] <- 0}  
	              else if(Z6[k] >= -0.67 & Z6[k] <= 0.67){Z.6[k] <- 1} 
                    else{Z.6[k] <- 2} 
                    }
       
                    Z = t(matrix(rbinom(n[i]*(p-6), 2, rep(q[-(1:6)], n[i])), ncol = n[i]))
                    X = cbind(Z.1, Z.2, Z.3, Z.4, Z.5, Z.6, Z)
                }
                for(t in 1:3){# loop over cor_y
                      if(t == 3){
                          mu = rep(0, n[i])
                      }
                      else if (t == 1){
                          mu = c(rep(5, floor(n[i]/3)), rep(1, floor(n[i]/3)), rep(5, n[i]-2*floor(n[i]/3)))
                      }
                      else{
                          mu = c(rep(4, floor(n[i]/3)), rep(3, floor(n[i]/3)), rep(2, n[i]-2*floor(n[i]/3)))
                      }
                       for(m in 1:3){
                            if(m ==3){ 
                                 Y = mu + 3 + X %*% beta.true + rnorm(n[i], 0, 1)
                            }
                            else if(m ==2){ 
                                 Y = mu + 3 + X %*% beta.true + X[,1]*X[,6] + X[,2]*X[,5] + 2* X[,3]^2 + rnorm(n[i], 0, 1)
                            }
                            else {
                                 Y = ((X[,1] + X[,2] + X[,3]) %% 2) * (X %*% beta.true) + rnorm(n[i], 0, 1)
                            }
                           idx = sample(7:p, h[j], replace = F)
                           for(c in 1:p){
                                if(c %in% idx){
                                id_c  = sample(1:n[i], floor(n[i]/2), rep = F) 
                                X[id_c, c] = X[id_c, c] + 1
                                X[-id_c, c] = X[-id_c, c] + 2                   
                                }
                           }
                           X = X%%3
                          ###single marker regression##
                          pValues=c(); beta.hat1 = c(); 
                          for(a in 1:p){
                                fm = lm(Y ~ X[,a])
						#print(c(s,i,j,l,t,m,a))
					  if(length(unique(X[,a])) == 1){
                                	pValues = c(pValues, 1)
                                	beta.hat1 = c(beta.hat1, 0)
					  }
					  if(length(unique(X[,a])) > 1){
                                	pValues = c(pValues, summary(fm)$coef[2,4])
                                	beta.hat1 = c(beta.hat1, summary(fm)$coef[2,1])
					  }							
                          }
                          index = order(pValues)
                          p.sort = sort(pValues)
                          if(length(which(p.sort < p.cutoff)) == 0 ){
                               Positive = 0
                          }
                          else{
                               Positive = index[which(p.sort < p.cutoff)]
                          }
                          fp1[s,count] = sum(Positive > 6)
                          tp1[s,count] = sum(Positive <= 6)
                          fdr1[s, count] = sum(Positive > 6)/length(Positive)
                          fpr1[s, count] = sum(Positive > 6)/(p - 6)
                          tpr1[s, count] = sum(Positive <= 6)/6
                          mse1[s,count] = var(beta.hat1) + sum((beta.hat1 - beta.true)^2)
                          marker1 = Positive
                          ##### LASSO with glmnet ######
	                   cv.fit = cv.glmnet(X, Y, alpha = 1)
	                   fit.lasso = glmnet(X, Y, family = "gaussian", alpha = 1, lambda = cv.fit$lambda.1se)
                         beta.hat2 = coef(fit.lasso)[-1]
                         fp2[s,count] = sum(beta.hat2[7:p] != 0)
                         tp2[s,count] = sum(beta.hat2[1:6] != 0)
                         s2[s,count] = sum(beta.hat2 != 0)
                         fdr2[s,count] = sum(beta.hat2[7:p] != 0) / max(s2[s,count],1)
                         tpr2[s,count] = sum(beta.hat2[1:6] != 0) / 6  
                         fpr2[s,count] = sum(beta.hat2[7:p] != 0) / (p - 6)  
                         mse2[s,count] = var(beta.hat2) + sum((beta.hat2 - beta.true)^2)
                         id.l = order(abs(beta.hat2), decreasing = TRUE)
                         marker2 = id.l[1:s2[s,count]]  
                         ##########BayesA########################################
                         ETA = list(MRK = list(X = X, model = 'BL'))
                         ETA$MRK$model = 'BL'
                         fmBA = BGLR(y = Y, ETA = ETA, nIter = nIter, burnIn = burnIn, saveAt = 'BL_')
                         beta.hat3 = fmBA$ETA[[1]]$b
                         fp3[s,count] = sum(abs(beta.hat3)[7:p] > beta3.cutoff)
                         tp3[s,count] = sum(abs(beta.hat3)[1:6] > beta3.cutoff)
                         s3[s,count] = sum(abs(beta.hat3) > beta.cutoff)
                         fdr3[s,count] = sum(abs(beta.hat3)[7:p] > beta.cutoff) / max(s3[s,count],1)
                         tpr3[s,count] = sum(abs(beta.hat3)[1:6] > beta.cutoff) / 6  
                         fpr3[s,count] = sum(abs(beta.hat3)[7:p] > beta.cutoff) / (p - 6)
                         mse3[s,count] = var(beta.hat3) + sum((beta.hat3 - beta.true)^2)
                         id.b = order(abs(beta.hat3), decreasing = TRUE)
                         marker3 = id.b[1:s3[s,count]]
                         nmv1[s,count] = sum(marker1 %in% marker2)
                         nmv2[s,count] = sum(marker1 %in% marker3)
                         nmv3[s,count] = sum(marker2 %in% marker3)
                         count = count + 1
				 #M = rbind(M, c(i,l,t,j,m))
                  }
              }
          }
       }
    }
}
#write.csv(M, file = "design.csv")
result1 = rbind(fdr1, fdr2, fdr3)
result2 = rbind(fpr1, fpr2, fpr3)
result3 = rbind(tpr1, tpr2, tpr3)
result4 = rbind(mse1, mse2, mse3)
result5 = rbind(nmv1, nmv2, nmv3)
result6 = rbind(fp1, fp2, fp3)
result7 = rbind(tp1, tp2, tp3)
save(result1, file = "fdr.RData")
save(result2, file = "fpr.RData")
save(result3, file = "tpr.RData")
save(result4, file = "mse.RData")
save(result5, file = "nmv.RData")
save(result6, file = "fp.RData")
save(result7, file = "tp.RData")