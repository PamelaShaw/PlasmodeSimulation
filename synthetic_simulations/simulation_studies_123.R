# Shaw, et al. A cautionary note for plasmode simulation studies in the setting of 
# causal inference. 

# Code for simulation studies 1, 2, 1a, 2a, 3
# Requires packages MatchIt, tmle, and their dependencies.

#------------------------------
# Function getEstimates
# Evaluate unadj, full PS Matching, IPTW, and TMLE
# estimates of the ATE based on correctly
# specified models for the PS and the outcome
# arguments: Y - outcome
# 			 A - binary treatment
# 			 W - baseline covariates
#	    family - regression family for outcome 
# returns: estimates from each approach
#------------------------------
getEstimates <- function(Y, A, W, Qform =  "Y ~ A + W1 + W2 + W3", gform = "A~W1 + W2 + W3", nfolds = 10, family = "gaussian"){
	require(MatchIt)
	require(tmle)
	m.unadj <- glm(Y ~ A, family = family)
	mu1.unadj <- predict(m.unadj, newdata = data.frame(A=1), type  = "response")
	mu0.unadj <- predict(m.unadj, newdata = data.frame(A = 0), type = "response")
	est.unadj <- mu1.unadj - mu0.unadj
#	matching
	if (length(Y) > 5000){
		method = "quick"
	} else {
		method = "full"
	}
 	matched <- suppressWarnings(MatchIt::matchit(formula = as.formula("A~W1 + W2 + W3"), 
                                                 data = data.frame(Y, A, W),
                                                 method = method, 
                                                 estimand = "ATE",
                                                 replace = TRUE))
 	mdat <- MatchIt::match.data(matched)
 #	linear regression gives identical ATE estimate as logistic regression when A is only predictor
   m.match <- glm(Y ~ A, data = mdat, weights = mdat$weights, family = family)
   mu1.match <- predict(m.match, newdata = data.frame(A = 1), type = "response")
   mu0.match <- predict(m.match, newdata = data.frame(A = 0), type = "response")

   est.match <- mu1.match - mu0.match
  
  
   
    ps <- predict(glm(gform, data = data.frame(A, W), family = "binomial"), type = "response")
    psSum <- summary(ps)
	wt.stab <- mean(A) * A/ps + mean(1-A) * (1-A) / (1-ps)
	wt.stab <- bound(wt.stab, c(0, 10^6))
	m.iptw.unbd <- glm(Y ~ A, weights = wt.stab, family = family)
	mu1.iptw.unbd <-predict(m.iptw.unbd, newdata = data.frame (A = 1), type = "response")
	mu0.iptw.unbd <-predict(m.iptw.unbd, newdata = data.frame (A = 0), type = "response")
	est.iptw.unbd <- mu1.iptw.unbd - mu0.iptw.unbd
	 

	
	wt.stab <- bound(wt.stab, c(0, sqrt(length(A)) * log(length(A)) / 5))
	m.iptw<- glm(Y ~ A, weights = wt.stab, family = family)
	mu1.iptw <-predict(m.iptw, newdata = data.frame (A = 1), type = "response")
	mu0.iptw  <-predict(m.iptw, newdata = data.frame (A = 0), type = "response")
	est.iptw  <- mu1.iptw - mu0.iptw
	 	  
	res.tmle <- tmle(Y, A, W, gform = gform ,Qform = Qform, family = family, evalATT = FALSE)
	res.tmle.unbd <- tmle(Y, A, W, gform = gform, Qform = Qform, family = family, gbound = 10^-6, evalATT = FALSE)
	
	m.glm <- glm(Qform, data = data.frame(Y, A, W), family = family)
	EY0.glm <- mean(predict(m.glm, newdata = data.frame(A=0, W), type ="response"))
	EY1.glm <- mean(predict(m.glm, newdata = data.frame(A=1, W), type ="response"))
	est.glm <- EY1.glm - EY0.glm
	
	m.glmPS <- glm(Y~ A + ps, data = data.frame(Y, A, W), family = family)
	EY0.glmPS <- mean(predict(m.glmPS, newdata = data.frame(A=0, ps), type ="response"))
	EY1.glmPS <- mean(predict(m.glmPS, newdata = data.frame(A=1, ps), type ="response"))
	est.glmPS <- EY1.glmPS - EY0.glmPS
	
	
	if(family == "binomial"){
		rr.unadj <- mu1.unadj/mu0.unadj
		rr.match <- mu1.match / mu0.match
	  	rr.iptw.unbd <- mu1.iptw.unbd / mu0.iptw.unbd
	  	rr.iptw  <- mu1.iptw / mu0.iptw
		rr.glm <- EY1.glm / EY0.glm
		rr.glmPS <- EY1.glmPS / EY0.glmPS
		
		RR = c(unadj = rr.unadj, match = rr.match, iptw = rr.iptw, iptw.unbd = rr.iptw.unbd, 
						tmle = res.tmle$estimates$RR$psi, tmle.unbd = res.tmle.unbd$estimates$RR$psi,
						rr.glm = rr.glm, rr.glmPS = rr.glmPS)
	} else {
		RR <- NULL
	}

	return(list(ATE = c(unadj = est.unadj, match = est.match, iptw = est.iptw, iptw.unbd = est.iptw.unbd,
						tmle = res.tmle$estimates$ATE$psi, tmle.unbd = res.tmle.unbd$estimates$ATE$psi,
						est.glm = est.glm, est.glmPS = est.glmPS), 
					RR = RR, 				
					psSum = psSum))
}

#------------------------------
# Demonstration that bootstrapping without modeling treatment can
# skew analysis results under correct model specification when
# there are near positivity violations.
# Two different DGPs
# 1. Continuous outcome  homogeneous treatment effect
# 2. Binary outcome homogeneous treatment effect
# Evaluate the ATE (risk difference when outcome is binary)

# For each DGP 
# Step 1. Define a data generating procedure, O = (W1, W2, W3, A, Y) 
# Step 2. Draw a sample with size n
# Step 3. 1000 bootstrap iterations 
#		a.One time only: Fit models for Y and A to the data using correctly specified parametric models
#		  because this investigation is about demonstrating the impact of retaining A vs. 
#		  randomly drawing A, not about bias do to unwarranted modeling assumptions.
#		REPEAT b-d 1000 TIMES
# 		b. Draw a sample of size n, O = (W,A), with replacement from the data generated in Step 2 
#		c. construct two datasets: 
#				d1 = (W,A,Y) based on original data (W,A), Y generated using Step 2 outcome model
#			 	d2 = (W, A',Y') based on original data (W), A', Y' generated using Step 2 PS and outcome model
# 		d. analyze both datasets to evaluate the ATE using unadjusted estimator, full PS matching, 
#			 IPTW with stabilized weights, TMLE, outcome regression, outcome regressing on the PS
# 			 all (except unadj) using correctly specified models
#		e. Evaluate bias, var, MSE of the 1000 bootstrap iterations based on d1 and d2
# Step 4. Compare estimators
#------------------------------

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


#------------------------------
# Step 1. Data Generating Processes
# Simulation study 1 (continuous),
# Simulation study 2 (binary)
# Simulation study 1a (continuous, treatment randomized)
# Simulation study 2a (binary, treatment randomized)
# Simulation study 3 (rare binary outcome)
#------------------------------
# Sim1
genSim.cont <- function(n, mult = 1.2, betaA = 2){
	W1 <- rnorm(n)
	W2 <- rnorm(n, mean = .5)
	W3 <- rbinom(n, 1, 0.4)
	A <- rbinom(n, 1, plogis(mult*(-.4 + .8*W1 + .01* W2 + .9*W3)))
	Y <- 10 + betaA*A + 1*W1 + 7*W2 - .02*W3 + rnorm(n)
	return(data.frame(Y, A, W1, W2, W3))
}
# Sim2
genSim.bin <- function(n, mult = 1.2, betaA = 1.1){
	W1 <- rnorm(n)
	W2 <- rnorm(n, mean = .5)
	W3 <- rbinom(n, 1, 0.4)
	A <- rbinom(n, 1, plogis(mult*(-.4 + .8*W1 + .01* W2 + .9*W3)))
	drs <- -1.8 + 0.24*W1 + .08*W2 + 0.8*W3 
	Y <-  rbinom(n, 1, plogis(betaA*A + drs))
	return(data.frame(Y, A, W1, W2, W3))
}
# Sim1a
genSim.cont.RCT <- function(n, mult = 1.2, betaA = 2){
	W1 <- rnorm(n)
	W2 <- rnorm(n, mean = .5)
	W3 <- rbinom(n, 1, 0.4)
	A <- rbinom(n, 1, .5)
	Y <- 10 + betaA*A + 1*W1 + 7*W2 - .02*W3 + rnorm(n)
	return(data.frame(Y, A, W1, W2, W3))
}
# Sim2a
genSim.bin.RCT <- function(n, mult = 1.2, betaA = 1.1){
	W1 <- rnorm(n)
	W2 <- rnorm(n, mean = .5)
	W3 <- rbinom(n, 1, 0.4)
	A <- rbinom(n, 1, .5)
	drs <- -1.8 + 0.24*W1 + .08*W2 + 0.8*W3 
	 Y <-  rbinom(n, 1, plogis(betaA*A + drs))
     return(data.frame(Y, A, W1, W2, W3))
}

# Sim3
genSim.bin.rare <- function(n, mult = .9, betaA = -2){
	W1 <- rnorm(n)
	W2 <- rnorm(n, mean = .5)
	W3 <- rbinom(n, 1, 0.4)
	A <- rbinom(n, 1, plogis(mult*(-.8 + .8*W1 + .2* W2 + .9*W3)))	
	drs <- -4.9+ 0.4*W1 - 4*W2 - 3*W3 
	Y <-  rbinom(n, 1, plogis(betaA*A + drs))
	return(data.frame(Y, A, W1, W2, W3))
}



# Step 2. Run the simulation study
#  Draw a sample (d.orig) that is the basis for a plasmode simulation study. 
#  Then run the plasmode sim two ways.
#    d1 is a simulated dataset based on (A,W) sampled from data
#    d2 is a simulated dataset based on (A.d2, W), where A.d2 is generated based on W.
runSim <- function(betaA = 2, B=1000, n = 100, nfolds = 20, family = "gaussian", 
	sim_fun, beta.Y, mult = 1.2){
    set.seed(1)
    d.orig <- sim_fun(n, mult = mult, betaA = betaA)
			
	# Step 3. Simulation study 1
	# 3a. fit correctly specified models to the data that will be used to generate A and Y
	m.A <- glm(A ~ ., data = d.orig[,-1], family = "binomial")
	pred.Y<- function(A, W){
			pred <-cbind(1,A,as.matrix(W))  %*%  beta.Y 
			if (family == "binomial"){
				pred <- plogis(pred)
			}
		return(pred)
	}
	
	 # simulated treatment effect comes from the fitted model
	 if(family == "gaussian"){
	 	psi0.sim <- beta.Y[2]
	 } else {
		psi0.sim <- mean(pred.Y(A=rep(1, nrow(d.orig)), W = d.orig[,-(1:2)])) - mean(pred.Y(A=rep(0, nrow(d.orig)), W = d.orig[,-(1:2)]))			
	}
	
	# 3b-d. run the study
	est.d1.ATE <- matrix(nrow = B, ncol = 8)
	colnames(est.d1.ATE) <- c("unadj", "Match", "IPTW", "IPTW.unbd", "TMLE", "TMLE.unbd", "glmCor", "glmPS")
	est.d2.ATE <-  est.d1.ATE
	psSum.d1 <- matrix(nrow = B, ncol = 6)
	colnames(psSum.d1) <- c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.")
	psSum.d2 <- psSum.d1
					
	for (b in 1:B){
		b.rows <- sample(1:n, n, replace = TRUE)
		# generate outcome using original A
		# generate random A for d2, A' = A.d2. Get model-based predicted outcome based on this new A.d2 and W
		A.d2 <- rbinom(n, 1, predict(m.A, newdata = d.orig[b.rows,], type = "response")) 
		pred.d1 <- pred.Y(A=d.orig$A[b.rows], d.orig[b.rows,-(1:2)])
		pred.d2 <- pred.Y(A=A.d2, d.orig[b.rows,-(1:2)])
		
		# Add the same random noise to the generated outcome in each analytic dataset
		if(family == "gaussian"){
			eps <- rnorm(n)
			Y.d1 <-  pred.d1 + eps
			Y.d2 <-  pred.d2 + eps
		}  else {
			U <- runif(n)
			Y.d1 <- pred.d1 >= U
			Y.d2 <- pred.d2 >= U
		}
		d1 <- getEstimates(Y = Y.d1, A = d.orig$A[b.rows], W = d.orig[b.rows,-(1:2)], nfolds = nfolds, family  = family)
		est.d1.ATE[b,] <- d1$ATE
	    psSum.d1[b,] <- d1$psSum
		d2 <- getEstimates(Y = Y.d2, A = A.d2, W = d.orig[b.rows,-(1:2)], nfolds = nfolds, family  = family)
		est.d2.ATE[b,] <- d2$ATE
		psSum.d2[b,] <- d2$psSum
	}
	
		# 3e. Evaluate results
	calcResults <- function(est, psi0){
		return(cbind(bias = colMeans(est - psi0),
						   medBias = apply(est-psi0, 2, median),
			   pctBias = colMeans(est - psi0) / psi0 * 100,
			   se = apply(est, 2, sd),
			   mse = colMeans((est - psi0)^2),
			   biasToSE = colMeans(est - psi0) / apply(est, 2, sd)
			   ))
	}	
	# Evaluate the true injected signal	
	res.d1.ATE <- calcResults(est.d1.ATE, psi0.sim)
	res.d2.ATE <- calcResults(est.d2.ATE, psi0.sim)		
		
	res.sim <- list(n = n, psSum.d1 = psSum.d1, psSum.d2 = psSum.d2,
				ATE= list(est.d1 = est.d1.ATE, est.d2 = est.d2.ATE,  
				res.d1=res.d1.ATE, res.d2=res.d2.ATE,psi0 = psi0.sim)) 
	return(res.sim)				
}

###

# Set number of bootstrap iterations and sample size
B <- 100000 
n <- 100


# Coefficient for generating data for Simulation study 1, 1a (continuous), 
# Simulation study 2, 2a (binary), Simulation study 3 (rare binary)
beta.Y.cont <- c(10, 2, 1, 7, .02)
beta.Y.bin <- c(-1.8, 1.1,  0.24, .08, .8)
beta.Y.bin.rare <- c(-4.9, -2, 0.4,  - 4, - 3)


# Run simulation studies	
res.sim1  <- runSim(B = B, n = n, family = "gaussian", sim_fun = genSim.cont, beta.Y = beta.Y.cont)
res.sim2  <- runSim(B = B, n = n, family = "binomial", sim_fun = genSim.bin, beta.Y = beta.Y.bin)
res.sim1a <- runSim(B = B, n = n, family = "gaussian", sim_fun = genSim.cont.RCT, beta.Y = beta.Y.cont)
res.sim2a <- runSim(B = B, n = n, family = "binomial", sim_fun = genSim.bin.RCT, beta.Y = beta.Y.bin)
res.sim3  <- runSim(B = B, n = n, family = "binomial", sim_fun = genSim.bin.rare, beta.Y = beta.Y.bin.rare, mult = 0.9)

 