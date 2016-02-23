#install.packages("bootstrap")
library(data.table)
library(bootstrap)
str(cholost)

d <- cholost$y # cholesterol decrease;
compliance <- cholost$z +0.0001 # compliance
c <- (compliance-mean(compliance))/sd(compliance) # standardized compliance
chol_data <- matrix(cbind(c, d), ncol=2)
n <- length(c)

set.seed(100)
#----------------------------------------------
# polynomial regression
#----------------------------------------------
c_first <- c[1]
sigma <- 22.0

table1_cp <- c()
for( i in 1:6){
	fit <- lm(d~poly(c, i, raw=TRUE))
	table1_cp[i] <- sum((fit$residual)^2) + 2*(sigma^2)*(i+1) - 80000
}

#table 1
cbind(2:7,table1_cp)
plot(c, d, xlim=c(-2,1.5), ylim=c(-50,120), main='figure1')
points(c, predict(lm(d~poly(c, 5, raw=TRUE))), col='blue')
first_fitted = predict( lm(d~poly(c, 5, raw=TRUE)), newdata=data.frame(c=c_first) )

#----------------------------------------------
# Bootstrap
#----------------------------------------------

B <- 4000
theta <- function(x, xdata){
	c_boot = xdata[x,1]	
	d_boot = xdata[x,2]
	
	table1_cp <- c()
	for( i in 1:6){
		fit <- lm(d_boot~poly(c_boot, i, raw=TRUE))
		table1_cp[i] <- sum((fit$residual)^2) + 2*(sigma^2)*(i+1) - 77000
	}
	Y=c()
	for( j in 1:length(d)){
		Y[j] = sum(d_boot == d[j])/sum(d == d[j])
	}
	index <- which(min(table1_cp) == table1_cp)
	fit <- lm(d_boot~poly(c_boot, index, raw=TRUE))
	pred <- predict( lm(d_boot~poly(c_boot, index, raw=TRUE)), newdata=data.frame(c_boot=c_first) )
	res <-  c(index, pred, Y)
	return(res)
}

results <- bootstrap(1:n, B, theta, chol_data)
res_table <- data.table(t(results$thetastar))
name = paste("V", 1:164, sep = "")
setnames(res_table,c('index', 'mean', name))

#table 1
100* table(res_table[,index]) /4000 # 2.87, 3.35, 2.85, 15.75, 43.67, 31.67
hist(res_table[,mean], xlab='bootstrap estimates for subject 1', breaks=30, main='figure3') #Figure 3
abline(v=first_fitted,col="red")
mean(res_table[,mean] < first_fitted) # 63.2% / 65.0%

#table 2
rbind(res_table[order(index) ,mean(mean), by = index][,V1]
,res_table[order(index) ,sd(mean), by = index][,V1])

#----------------------------------------------
# Accuracy of the smoothed bootstrap estimates
#----------------------------------------------

prop_est <- mean(res_table[,mean])
T <- res_table[,mean] - prop_est
prop_var <- 0
for( str in name){
	prop_var <- prop_var + (mean(res_table[,get(str)]*T))^2
}
prop_sd <- sqrt(prop_var)

#table 3
standard_ci <- c( first_fitted + qnorm(0.025)*sd(res_table[,mean]), first_fitted - qnorm(0.025)*sd(res_table[,mean]), first_fitted )
percentile_ci <- c( quantile(x = res_table[,mean], probs=c(0.025, 0.975)), mean( quantile(x = res_table[,mean], probs=c(0.025, 0.975)) ) )
prop_ci <- c( prop_est + qnorm(0.025)*prop_sd, prop_est - qnorm(0.025)*prop_sd, prop_est )

rbind( standard_ci, percentile_ci, prop_ci)



