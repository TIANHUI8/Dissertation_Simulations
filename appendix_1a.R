setwd("D:/Clemson/RESEARCH/RESEARCH/mypaper/simulation")
r1 = seq(0, 1, by = 0.01)
r2 = seq(0, 1, by = 0.01) + 0.01
rep = 1000
n = 100
p = 50
lenr = length(r1)
sum = matrix(-99, ncol = rep, nrow = lenr)
q = runif(p, min = 0, max = 1)
for (i in 1:lenr){
      for (j in 1:rep){
		X = t(matrix(rbinom(n*p, 2, rep(q, n)), ncol = n))
		Cor = cor(X)[1, 2:p]
		Cor[is.na(Cor)] = 0
            c1 = rep(r1[i], p - 1)
            c2 = rep(r2[i], p - 1)
		sum[i,j] =sum(sum(abs(Cor) <= c2 & abs(Cor) > c1) > 0)
		}
}
prob = apply(sum, 1, mean)
write.csv(prob, file = "appendix1_Prob5.csv")
#############read result ########################
prob1 = read.csv(file = "appendix1_Prob5.csv")[,-1]
prob2 = read.csv(file = "appendix1_Prob6.csv")[,-1]
prob3 = read.csv(file = "appendix1_Prob7.csv")[,-1]
prob4 = read.csv(file = "appendix1_Prob8.csv")[,-1]
prob1 = stepfun(r1[-1], prob1, f = 0)
prob2 = stepfun(r1[-1], prob2, f = 0)
prob3 = stepfun(r1[-1], prob3, f = 0)
prob4 = stepfun(r1[-1], prob4, f = 0)
par(mfrow = c(2,2))
plot(prob1, lty = 3, col = "red", xlab = "r", ylab = "probability", ylim = c(0, 1), lwd = 2, pch = 1, main = "n = 100, p = 50")
plot(prob2, lty = 3, col = "red", xlab = "r", ylab = "probability", ylim = c(0, 1), lwd = 2, pch = 1, main = "n = 100, p = 100")
plot(prob3, lty = 3, col = "red", xlab = "r", ylab = "probability", ylim = c(0, 1), lwd = 2, pch = 1, main = "n = 100, p = 500")
plot(prob4, lty = 3, col = "red", xlab = "r", ylab = "probability", ylim = c(0, 1), lwd = 2, pch = 1, main = "n = 100, p = 1000")
