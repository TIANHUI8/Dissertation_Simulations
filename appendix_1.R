setwd("D:/Clemson/RESEARCH/RESEARCH/mypaper/simulation")
r1 = 0.2
r2 = 0.6
rep = 1000
n = seq(10, 100, by = 10)
p = seq(50, 1000, by = 50)
lenn = length(n)
lenp = length(p)
sum = matrix(0, ncol = lenp, nrow = lenn)
for (j in 1:lenp){
      q = runif(p[j], min = 0, max = 1)
      for(i in 1:lenn){
		for (k in 1:rep){
		   	X = t(matrix(rbinom(n[i]*p[j], 2, rep(q, n[i])), ncol = n[i]))
			Cor = cor(X)[1, 2:p[j]]
			Cor[is.na(Cor)] = 0
                  c1 = rep(r1, p[j] - 1)
			c2 = rep(r2, p[j] - 1)
			sum[i,j] = sum[i,j] + sum(sum(abs(Cor) <= c2 & abs(Cor) > c1) > 0)
		}
	}
}
prob = sum / rep
write.csv(prob, file = "appendix1_probb.csv")
#############read result ########################
prob = read.csv(file = "appendix1_prob.csv")[,-1]
prob1 = read.csv(file = "appendix1_prob1.csv")[,-1]
prob2 = read.csv(file = "appendix1_prob2.csv")[,-1]
prob3 = read.csv(file = "appendix1_prob3.csv")[,-1]
prob4 = read.csv(file = "appendix1_prob4.csv")[,-1]
par(mfrow = c(2,2))
plot(n, prob1[,20], type = "o", lty = 3, col = "red", ylab = "probability", ylim = c(0, 1), lwd = 2, pch = 1, main = "(0,0.25], p = 1000")
plot(n, prob2[,20], type = "o", lty = 3, col = "red", ylab = "probability", ylim = c(0, 1), lwd = 2, pch = 1, main = "(0.25,0.5], p = 1000")
plot(n, prob3[,20], type = "o", lty = 3, col = "red", ylab = "probability", ylim = c(0, 1), lwd = 2, pch = 1, main = "(0.5,0.75], p = 1000")
plot(n, prob4[,20], type = "o", lty = 3, col = "red", ylab = "probability", ylim = c(0, 1), lwd = 2, pch = 1, main = "(0.75,1], p = 1000")
par(mfrow = c(2,2))
plot(n, prob[,1], type ="o", lty = 3, col = "red", ylab = "probability", ylim = c(0, 1), lwd = 2, pch = 1, main = "[0.2,0.6], p = 50")
plot(n, prob[,10], type = "o", lty = 3, col = "red", ylab = "probability", ylim = c(0, 1), lwd = 2, pch = 1, main = "[0.2,0.6], p = 500")
plot(n, prob[,15], type = "o", lty = 3, col = "red", ylab = "probability", ylim = c(0, 1), lwd = 2, pch = 1, main = "[0.2,0.6], p = 750")
plot(n, prob[,20], type = "o", lty = 3, col = "red", ylab = "probability", ylim = c(0, 1), lwd = 2, pch = 1, main = "[0.2,0.6], p = 1000")
