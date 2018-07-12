# Things to change
path = "/scratch2/tianhui/tianhui/GWAS/"		## path of results folders
len = 240				## number of folders
sims = 5				## number of simulations per folder
FDR1 = matrix(-99, nrow = len*sims, ncol = 243)
FDR2 = matrix(-99, nrow = len*sims, ncol = 243)
FDR3 = matrix(-99, nrow = len*sims, ncol = 243)
FPR1 = matrix(-99, nrow = len*sims, ncol = 243)
FPR2 = matrix(-99, nrow = len*sims, ncol = 243)
FPR3 = matrix(-99, nrow = len*sims, ncol = 243)
TPR1 = matrix(-99, nrow = len*sims, ncol = 243)
TPR2 = matrix(-99, nrow = len*sims, ncol = 243)
TPR3 = matrix(-99, nrow = len*sims, ncol = 243)
MSE1 = matrix(-99, nrow = len*sims, ncol = 243)
MSE2 = matrix(-99, nrow = len*sims, ncol = 243)
MSE3 = matrix(-99, nrow = len*sims, ncol = 243)
NMV1 = matrix(-99, nrow = len*sims, ncol = 243)
NMV2 = matrix(-99, nrow = len*sims, ncol = 243)
NMV3 = matrix(-99, nrow = len*sims, ncol = 243)
id.error = c()

for(i in 1:len){

	id = ((i-1)*sims + 1):(i*sims)
	res1 = tryCatch(file.exists(paste0(path,i,"/fdr.RData")), error = function(e) 0)  
	res2 = tryCatch(file.exists(paste0(path,i,"/fpr.RData")), error = function(e) 0)	
	res3 = tryCatch(file.exists(paste0(path,i,"/tpr.RData")), error = function(e) 0)
        res4 = tryCatch(file.exists(paste0(path,i,"/mse.RData")), error = function(e) 0)
        res5 = tryCatch(file.exists(paste0(path,i,"/nmv.RData")), error = function(e) 0)
      
	if(res1 == 0 | res2 == 0 | res3 == 0 | res4 == 0 |res5 == 0){
		id.error = c(id.error, id)
		next
	}
  	load(paste0(path,i,"/fdr.RData"))
	load(paste0(path,i,"/fpr.RData"))
	load(paste0(path,i,"/tpr.RData"))
        load(paste0(path,i,"/mse.RData"))
        load(paste0(path,i,"/nmv.RData"))

	## Things to change
	FDR1[id,] = result1[1:sims,]
	FDR2[id,] = result1[(sims+1):(2*sims),]
	FDR3[id,] = result1[(2*sims+1):(3*sims),]
        FPR1[id,] = result2[1:sims,]
        FPR2[id,] = result2[(sims+1):(2*sims),]
        FPR3[id,] = result2[(2*sims+1):(3*sims),]
        TPR1[id,] = result3[1:sims,]
        TPR2[id,] = result3[(sims+1):(2*sims),]
        TPR3[id,] = result3[(2*sims+1):(3*sims),]
        MSE1[id,] = result4[1:sims,]
        MSE2[id,] = result4[(sims+1):(2*sims),]
        MSE3[id,] = result4[(2*sims+1):(3*sims),]
        NMV1[id,] = result5[1:sims,]
        NMV2[id,] = result5[(sims+1):(2*sims),]
        NMV3[id,] = result5[(2*sims+1):(3*sims),] 
}

## Things to change
if(length(id.error) > 0){
	FDR1 = FDR1[-id.error,]
        FDR2 = FDR2[-id.error,]
        FDR3 = FDR3[-id.error,]
	FPR1 = FPR1[-id.error,]
        FPR2 = FPR2[-id.error,]
        FPR3 = FPR3[-id.error,]
        TPR1 = TPR1[-id.error,]
        TPR2 = TPR2[-id.error,] 
        TPR3 = TPR3[-id.error,] 
        MSE1 = MSE1[-id.error,] 
        MSE2 = MSE2[-id.error,]
        MSE3 = MSE3[-id.error,]
        NMV1 = NMV1[-id.error,] 
        NMV2 = NMV2[-id.error,]
        NMV3 = NMV3[-id.error,]

}

## Things to change
res1 = list("fdr" = FDR1, "fpr" = FPR1, "tpr" = TPR1, "mse" = MSE1, "nmv" = NMV1)
res2 = list("fdr" = FDR2, "fpr" = FPR2, "tpr" = TPR2, "mse" = MSE2, "nmv" = NMV2)
res3 = list("fdr" = FDR3, "fpr" = FPR3, "tpr" = TPR3, "mse" = MSE3, "nmv" = NMV3)


save(res1, file = "results1.merged.RData")
save(res2, file = "results2.merged.RData")
save(res3, file = "results3.merged.RData")

save(id.error, file = "id.error.RData")
