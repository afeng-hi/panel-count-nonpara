library(splines2)


len1 = INN1 + nknots1 - 1
len2 = INN2 + nknots2 - 1
total_len = len1 + len2 + length(alpha) 


id.alpha = (len1 + len2 + 1) : total_len


id = which(!is.na(theta.est[, id.alpha[1]]) & abs(theta.est[, id.alpha[1]] - alpha[1]) < 3 &
             !is.na(theta.est[, id.alpha[2]]) & abs(theta.est[, id.alpha[2]] - alpha[2]) < 3)
id =1:200
cat("============== alpha ==============\n")
cat("sim:", length(id), "\n\n")


for(k in 1:length(alpha)){
  cat(">>> alpha_", k, " (true: ", alpha[k], ") <<<\n", sep="")
  
  bias = mean(theta.est[id, id.alpha[k]]) - alpha[k]
  sd_val = sd(theta.est[id, id.alpha[k]])
  ese_val = mean(bots.sd[id, id.alpha[k]], na.rm=TRUE)

  JL = (theta.est[id, id.alpha[k]] - 1.96 * bots.sd[id, id.alpha[k]]) < alpha[k]
  JR = (theta.est[id, id.alpha[k]] + 1.96 * bots.sd[id, id.alpha[k]]) > alpha[k]
  cp_val = mean(JL * JR, na.rm=TRUE)
  
  cat("Bias:", bias, "\n")
  cat("SD  :", sd_val, "\n")
  cat("ESE :", ese_val, "\n")
  cat("CP   (95%):", cp_val, "\n\n")
}
cat("========================================================\n")



id.xi1 = 1 : len1
id.xi2 = (len1 + 1) : (len1 + len2)

thetahat_lam = apply(theta.est[id, id.xi1], 2, mean, na.rm=TRUE)
thetahat_x   = apply(theta.est[id, id.xi2], 2, mean, na.rm=TRUE)




len_pts = 3000


x_lam = seq(0, tau, length.out = len_pts)
BB1 = bSpline(x_lam, df=INN1, knots=knots1, degree=(INN1-1), intercept=FALSE)

y_lam = BB1 %*% thetahat_lam  
true_lam = x_lam             


all_y_lam = BB1 %*% t(theta.est[id, id.xi1]) 

lower_lam = apply(all_y_lam, 1, quantile, probs = 0.025, na.rm = TRUE)
upper_lam = apply(all_y_lam, 1, quantile, probs = 0.975, na.rm = TRUE)



x_x = seq(0, 1, length.out = len_pts)
BB2 = bSpline(x_x, df=INN2, knots=knots2, degree=(INN2-1), intercept=FALSE)

y_x = BB2 %*% thetahat_x      
true_x = 2 * x_x^2            


all_y_x = BB2 %*% t(theta.est[id, id.xi2]) 
lower_x = apply(all_y_x, 1, quantile, probs = 0.025, na.rm = TRUE)
upper_x = apply(all_y_x, 1, quantile, probs = 0.975, na.rm = TRUE)






# ===================== Lambda_0(s) 画图 =====================
plot(x_lam, true_lam, type="l", col="black", lty=1, lwd=2, 
     ylim=range(c(lower_lam, upper_lam, true_lam, y_lam), na.rm=TRUE),
     ylab=expression(Lambda), xlab="s", main="n=400 and censoring rate=40%")


lines(x_lam, y_lam, type="l", col="red", lty=2, lwd=2)


lines(x_lam, lower_lam, type="l", col="blue", lty=4, lwd=1.5)
lines(x_lam, upper_lam, type="l", col="blue", lty=4, lwd=1.5)


legend("topleft", lty=c(1, 2, 4), col=c("black", "red", "blue"), 
       lwd=c(2, 2, 1.5), legend=c("True Function", "Mean of the Estimates", "Pointwise Percentiles"))


# ===================== beta(x) 画图 =====================
plot(x_x, true_x, type="l", col="black", lty=1, lwd=2, 
     ylim=range(c(lower_x, upper_x, true_x, y_x), na.rm=TRUE),
     ylab=expression(beta), xlab="x", main="n=400 and censoring rate=40%")


lines(x_x, y_x, type="l", col="red", lty=2, lwd=2)


lines(x_x, lower_x, type="l", col="blue", lty=4, lwd=1.5)
lines(x_x, upper_x, type="l", col="blue", lty=4, lwd=1.5)


legend("topleft", lty=c(1, 2, 4), col=c("black", "red", "blue"), 
       lwd=c(2, 2, 1.5), legend=c("True Function", "Mean of the Estimates", "Pointwise Percentiles"))