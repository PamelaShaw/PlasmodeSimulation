# Shaw, et al. A cautionary note for plasmode simulation studies in the setting of 
# causal inference. 
# 
# Code for simulation studies 4a and 4b
# Estimating MSM parameter

getIPTW_MSM <- function(Y, A, W, msm = "Y ~ A + W1 + W2 + W3",  gform, family = "binomial"){

    ps <- predict(glm(gform, data = data.frame(A, W), family = "binomial"), type = "response")
    psSum <- summary(ps)
	wt.stab <- mean(A) * A/ps + mean(1-A) * (1-A) / (1-ps)
	wt.stab <- bound(wt.stab, c(0, sqrt(length(A)) * log(length(A)) / 5))
	if(family == "binomial"){
		family = "quasibinomial"
	}
	est.iptw.msm <- coef(glm(msm, data = data.frame(Y, A, W, weights = wt.stab), weights = weights, family = family))
		 
	 	return(list(msm = est.iptw.msm, 
					psSum = psSum))
}


#------------------------------
# function bound: Utility function for
# truncating values to lie between
# the specified min and max bounds
#------------------------------
bound <- function(x, bds){
	x[x<min(bds)] <- min(bds)
	x[x>max(bds)] <- max(bds)
	return(x)
} 


# Data generating process
genSim.bin.msm <- function(n, mult = 1, betaA = 1.1, beta.W = c(-2.5, 0.24, 0.08, 0.8, -0.3, -0.6)){
	W1 <- rnorm(n)
	W2 <- rnorm(n, mean = .5)
	W3 <- rbinom(n, 1, 0.4)
	W4 <- rnorm(n) + W3
	W5 <- as.integer(W1 > .2)
	
	A <- rbinom(n, 1,  plogis(mult*(-.4 + .8*W1 + .01* W2 + .9*W3)))
	
	logitY <- cbind(1,A,W1, W2, W3, W4, W5)  %*%  c(beta.Y[1], betaA, beta.W[-1]) # not a typo
	Y <-  rbinom(n, 1, plogis(logitY))     # approx 14%
	 
	return(data.frame(Y, A, W1, W2, W3, W4, W5))
}



runSim.msm <- function(B=1000, n = 100, msm = "Y ~ A + W1 + W2 + W3", 
  family = "binomial", beta.Y, gform = "A ~ W1 + W2 + W3"){
    set.seed(1)
	d.orig <- genSim.bin.msm(n)
	# Step 3. Simulation study 1
	# 3a. fit correctly specified model to the data that will be used to generate A 
	m.A <- glm(gform, data = d.orig[,-1], family = "binomial")
	# model Y the same way every time.
	pred.Y<- function(A, W, family){
		pred <-cbind(1,A,as.matrix(W))  %*%  beta.Y 
		if(family == "binomial"){
			pred <- plogis(pred)
		}
		return(pred)
	}
	# the true parameter is the projection of the truth onto the msm in d.orig
	psi0.msm<- getIPTW_MSM(Y = d.orig$Y, A = d.orig$A, W = d.orig[,-(1:2)], msm = msm, gform = gform, family  = family)$msm

	
	# 3b-d. run the study
	est.d1.msm <- est.d2.msm <- matrix(NA, nrow = B, ncol = 7)
	colnames(est.d1.msm) <- colnames(est.d2.msm) <- c("int", "A", "W1", "W2", "W3", "W4", "W5")
	psSum.d1 <- matrix(nrow = B, ncol = 6)
	colnames(psSum.d1) <- c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.")
	psSum.d2 <- psSum.d1
	set.seed(10) 
	for (b in 1:B){
		b.rows <- sample(1:n, n, replace = TRUE)
		# for dataset d1 get predicted probability outcome = 1 based on original A, W
		# for dataset d2 get predicted value, A.d2, conditional on W, then
		# predicted probability outcome = 1 using A.d2, W
		A.d2 <- rbinom(n, 1, predict(m.A, newdata = d.orig[b.rows,], type = "response")) 
		pred.d1 <- pred.Y(A=d.orig$A[b.rows], d.orig[b.rows,-(1:2)], family = family)
		pred.d2 <- pred.Y(A=A.d2, d.orig[b.rows,-(1:2)], family = family)

		# Avoid datasets in which there are no events in one of the arms, 
		# If that happens in either d1 or d2 re-draw the outcome from the Bernoulli distribution again  
		noEvents <- TRUE  
		while(noEvents){
			U <- runif(n)
			Y.d1 <- pred.d1 >= U
			Y.d2 <- pred.d2 >= U
			noEvents <- any(c(table(Y.d1, d.orig$A[b.rows]), table(Y.d2, A.d2)) == 0)
		}
		
		# Evaluate MSM parameter
		d1 <- getIPTW_MSM(Y = Y.d1, A = d.orig$A[b.rows], W = d.orig[b.rows,-(1:2)], msm = msm, gform = gform, family  = family)
		est.d1.msm[b,1:length(d1$msm)] <- d1$msm
		psSum.d1[b,] <-  d1$psSum
		
		d2 <- getIPTW_MSM(Y = Y.d2, A = A.d2, W = d.orig[b.rows,-(1:2)], msm = msm, gform = gform, family  = family)
		est.d2.msm[b,1:length(d2$msm)] <- d2$msm
		psSum.d2[b,] <- d2$psSum			
	}
	return(list(psi0.msm = psi0.msm, est.d1.msm=est.d1.msm,  est.d2.msm=est.d2.msm, 
						psSum.d1 = psSum.d1, psSum.d2 = psSum.d2))
}
		
B <- 10000
n <- 100	

beta.Y <- c(-2.5, 1.1, 0.24, 0.08, 0.8, -0.3, -0.6)
beta.W <- beta.Y[-2]

res.sim4a <- runSim.msm(B = B, n = n, msm = "Y ~ A + W1 + W2 + W3 + W4 + W5", family = "binomial", beta.Y = beta.Y)
res.sim4b <- runSim.msm(B = B, n = n, msm = "Y ~ A + W1 + W2 + W3", family = "binomial", beta.Y = beta.Y)
	 
