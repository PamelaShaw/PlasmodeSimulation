# utility functions for plasmode paper simulation
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
# Function getEstimates
# Evalute unadj, full PS Matching, IPTW, and TMLE
# estimates of the ATE based on correctly
# specified models for the PS and the outcome
# arguments: Y - outcome
# 			 A - binary treatment
# 			 W - baseline covariates
#	    family - regression family for outcome 
# returns: estimates from each approach
#------------------------------
getEstimates <- function(Y, A, W, Qform =  "Y ~ A + W1 + W2 + W3", 
                         gform = "A~W1 + W2 + W3", nfolds = 10, 
                         family = "gaussian"){
  require(MatchIt)
  require(tmle)
  #require(hal9001)
  psStrat <- function(ps, nstrata = 10, param = "ATE"){
    quantiles <- seq(0, 1, length = nstrata) 
    ps.strata <- cut(ps, breaks=unique(quantile(ps, p=quantiles)), include.lowest=TRUE, labels=FALSE)
    psi.strat <- by(cbind(Y,A), ps.strata,function(x){coef(lm(x[,1] ~ x[,2], data = as.data.frame(x)))[2]})
    if (param == "ATE") {
      wt <- table(ps.strata)/length(ps)
    } else {
      # ATT weight is proportion treated in each stratum
      wt <- by(cbind(A), ps.strata, sum)/sum(A)
    }
    return(sum(wt * psi.strat))
  }
  unadj_fit <- lm(Y ~ A)
  est.unadj <- coef(unadj_fit)[2]	
  ey0_unadj <- mean(predict(unadj_fit, newdata = data.frame(A = 0, W), type = "response"))
  ey1_unadj <- mean(predict(unadj_fit, newdata = data.frame(A = 1, W), type = "response"))
  #	matching
  if (length(Y) > 5000){
    method = "quick"
  } else {
    method = "full"
  }
  # matching ATE, EY0, EY1
  matched <- suppressWarnings(MatchIt::matchit(formula = as.formula(gform), 
                                               data = data.frame(Y, A, W),
                                               method = method, 
                                               estimand = "ATE",
                                               replace = TRUE))
  mdat <- MatchIt::match.data(matched)
  #	linear regression gives identical ATE estimate as logistic regression when A is only predictor
  df_0 <- data.frame(A = 0, W)
  df_1 <- data.frame(A = 1, W)
  match_fit <- lm(Y ~ A, data = mdat, weights = mdat$weights)
  ey0_match <- mean(predict(match_fit, newdata = df_0, type = "response"))
  ey1_match <- mean(predict(match_fit, newdata = df_1, type = "response"))
  est.match <- coef(match_fit)[2]    
  # IPTW
  ps <- predict(glm(gform, data = data.frame(A, W), family = "binomial"), 
                type = "response")
  psSum <- summary(ps)
  wt.stab.unbd <- mean(A) * A/ps + mean(1-A) * (1-A) / (1-ps)
  wt.stab.unbd <- bound(wt.stab.unbd, c(0, 10^6))
  iptw_fit_unbd <- lm(Y ~ A, weights = wt.stab.unbd)
  est.iptw.unbd <- coef(iptw_fit_unbd)[2] 
  wt.stab <- bound(wt.stab.unbd, c(0, sqrt(length(A)) * log(length(A)) / 5))
  iptw_fit <- lm(Y ~ A, weights = wt.stab)
  est.iptw <- coef(iptw_fit)[2]
  ey0_iptw_unbd <- mean(predict(iptw_fit_unbd, newdata = df_0, type = "response"))
  ey1_iptw_unbd <- mean(predict(iptw_fit_unbd, newdata = df_1, type = "response"))
  ey0_iptw <- mean(predict(iptw_fit, newdata = df_0, type = "response"))
  ey1_iptw <- mean(predict(iptw_fit, newdata = df_1, type = "response"))
  # TMLE
  res.tmle <- tmle(Y, A, W, gform = gform, Qform = Qform, family = family)
  res.tmle.unbd <- tmle(Y, A, W, gform = gform, Qform = Qform, family = family, 
                        gbound = 10^-6)
  ey0_tmle <- res.tmle$estimates$EY0$psi
  ey1_tmle <- res.tmle$estimates$EY1$psi
  ey0_tmle_unbd <- res.tmle.unbd$estimates$EY0$psi
  ey1_tmle_unbd <- res.tmle.unbd$estimates$EY1$psi
  # GLMs
  m.glm <- glm(Qform, data = data.frame(Y, A, W), family = family)
  ey0_glm <- mean(predict(m.glm, newdata = data.frame(A=0, W), type ="response"))
  ey1_glm <- mean(predict(m.glm, newdata = data.frame(A=1, W), type ="response"))
  est.glm <- ey1_glm - ey0_glm
  cate_glm <- coefficients(m.glm)[2]
  m.glmPS <- glm(Y~ A + ps, data = data.frame(Y, A, W), family = family)
  ey0_ps <- mean(predict(m.glmPS, newdata = data.frame(A=0, ps), type ="response"))
  ey1_ps <- mean(predict(m.glmPS, newdata = data.frame(A=1, ps), type ="response"))
  est.glmPS <- ey1_ps - ey0_ps
  est.PSstrat5 <- psStrat(ps, nstrata = 5)
  est.PSstrat10 <- psStrat(ps, nstrata = 10)
  est.PSstrat20 <- psStrat(ps, nstrata = 20)
  
  EY0 <- c(unadj = ey0_unadj, match = ey0_match, iptw = ey0_iptw,
           iptw.unbd = ey0_iptw_unbd, tmle = ey0_tmle,
           tmle.unbd = ey0_tmle_unbd, est.glm = ey0_glm,
           est.glmPS = ey0_ps)
  EY1 <- c(unadj = ey1_unadj, match = ey1_match, iptw = ey1_iptw,
           iptw.unbd = ey1_iptw_unbd, tmle = ey1_tmle,
           tmle.unbd = ey1_tmle_unbd, est.glm = ey1_glm,
           est.glmPS = ey1_ps)
  
  return(list(ATE = c(unadj = est.unadj, match = est.match, iptw = est.iptw, 
                      iptw.unbd = est.iptw.unbd,
                      tmle = res.tmle$estimates$ATE$psi,
                      tmle.unbd = res.tmle.unbd$estimates$ATE$psi,
                      est.glm = est.glm, est.glmPS = est.glmPS, 
                      est.PSstrat5 = est.PSstrat5, est.PSstrat10 = est.PSstrat10, 
                      est.PSstrat20 = est.PSstrat20), 
              psSum = psSum,
              EY0 = EY0,
              EY1 = EY1,
              RR = EY1 / EY0,
              CATE = c(glm = cate_glm)))
}

#------------------------------
# Step 1. Function genSim.cont and genSim.bin
# 	generate treatment as outcome as a function of 
#   two normally distributed and 1 binary covariate
# 	the PS in a sample of size 10^6 is between 0.014 and 0.990
#   so positivity issues when n < 2000 or so (based on our 
#   truncation bound of 5/sqrt(n)/ln(n) = 0.015 at n=2000)
#   the true additive treatment effect = 2.
# when mult =1 ps between 0.05 and 0.95
# when mult = 1.2, n = 10^5 
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.01333 0.31158 0.49583 0.49723 0.68120 0.98864
# when mult = 1.5 ps is  
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.005751 0.271626 0.497086 0.496932 0.724140 0.996596 
#------------------------------
genSim.cont <- function(n, mult = 1.2, betaA = 2){
  W1 <- rnorm(n)
  W2 <- rnorm(n, mean = .5)
  W3 <- rbinom(n, 1, 0.4)
  A <- rbinom(n, 1, plogis(mult*(-.4 + .8*W1 + .01* W2 + .9*W3)))
  Y <- 10 + betaA*A + 1*W1 + 7*W2 - .02*W3 + rnorm(n)
  return(data.frame(Y, A, W1, W2, W3))
}

genSim.bin <- function(n, mult = 1.2, betaA = 1.1){
  W1 <- rnorm(n)
  W2 <- rnorm(n, mean = .5)
  W3 <- rbinom(n, 1, 0.4)
  A <- rbinom(n, 1, plogis(mult*(-.4 + .8*W1 + .01* W2 + .9*W3)))
  drs <- -1.8 + 0.24*W1 + .08*W2 + 0.8*W3 
  Y <-  rbinom(n, 1, plogis(betaA*A + drs))
  #	mu0 <- mean(plogis(drs))
  #   mu1 <- mean(plogis(betaA*A + drs))
  return(data.frame(Y, A, W1, W2, W3))
}

# @param original_covariates a data.frame with the original covariates (W)
# @param n the desired sample size
# @param ps_model a regression model for the propensity of treatment given W
# @param outcome_model a regression model for the outcome given treatment and W
# @return a list of two plasmode datasets: one with (A, W) sampled together 
#         and Y generated from there; and one with W sampled, then A | W, 
#         then Y | A,W
gen_plasmode_dataset <- function(original_covariates = NULL, 
                                 n = nrow(original_covariates), 
                                 original_tx = NULL, ps_model = NULL, 
                                 outcome_model = NULL) {
  # get the bootstrap rows
  b_rows <- sample(1:nrow(original_covariates), n, replace = TRUE)
  new_covariates <- original_covariates[b_rows, ]
  # generate new A
  a <- rbinom(n, size = 1, prob = predict(ps_model, newdata = new_covariates, 
                                          type = "response"))
  # computed E(Y | A,W) based on A, W sampled together
  pred_y_1 <- pred.Y(A = original_tx[b_rows], W = new_covariates, 
                     outcome_model = outcome_model)
  # compute E(Y | A,W) based on A | W, W 
  pred_y_2 <- pred.Y(A = a, W = new_covariates, outcome_model = outcome_model)
  # generate binary Y
  U <- runif(n)
  y_1 <- as.numeric(pred_y_1 >= U)
  y_2 <- as.numeric(pred_y_2 >= U)
  # return datasets
  return(list("df_aw" = data.frame(Y = y_1, A = original_tx[b_rows], new_covariates),
              "df_w" = data.frame(Y = y_2, A = a, new_covariates),
              "pred_y_aw" = pred_y_1, "pred_y_w" = pred_y_2))
}

# pred.Y <- function(A, W, outcome_model, beta.Y){
    # pred <- cbind(1,A,as.matrix(W))  %*%  beta.Y 
    # model_matrix <- model.matrix(outcome_model$formula, 
    #                              data = data.frame(SH_HOSP_1826days = 1,
    #                                                Y = 1, 
    #                                                numEpiType = A, A = 1, W))
    # pred <- model_matrix %*% beta.Y
    # #if (family == "binomial"){
    #   pred <- plogis(pred)
    # #}
pred.Y <- function(A, W, outcome_model) {
    pred <- predict(outcome_model, newdata = data.frame(A = A, W), type = "response")
    return(pred)
  }
  

# summarize results
# @param est the point estimates
# @param psi0 the estimand
calcResults <- function(est, psi0, na.rm = FALSE){
  return(cbind(bias = colMeans(est - psi0),
               medBias = apply(est-psi0, 2, median),
               pctBias = colMeans(est - psi0) / psi0 * 100,
               se = apply(est, 2, sd),
               mse = colMeans((est - psi0)^2),
               biasToSE = colMeans(est - psi0) / apply(est, 2, sd)
  ))
}	

# run one iteration of the simulation
# @param b the iteration number
# @param n the sample size
# @param original_covariates the original covariates
# @param original_tx the original treatment assignment
# @param ps_model the treatment propensity score model fit to the original data
# @param outcome_model the outcome model fit to the original data
# @param family the family
# @return the estimates fit to both datasets
one_sim <- function(b = 1, original_covariates = NULL, 
                    n = nrow(original_covariates),
                    original_tx = NULL, ps_model = NULL, outcome_model = NULL,
                    family = "binomial") {
  # get plasmode datasets
  datasets <- gen_plasmode_dataset(original_covariates = original_covariates, 
                                   n = n, original_tx = original_tx, 
                                   ps_model = ps_model, 
                                   outcome_model = outcome_model)
  # get correlation between A, W in each dataset
  cor_aw <- cor(datasets$df_aw$A, datasets$df_aw[, -c(1:2)])
  cor_w <- cor(datasets$df_w$A, datasets$df_w[, -c(1:2)])
  # set up
  nice_qform <- paste0(format(formula(outcome_model)), collapse = "")
  nice_gform <- paste0(format(formula(ps_model)), collapse = "")
  # fit estimators
  ests_aw <- getEstimates(Y = datasets$df_aw$Y, A = datasets$df_aw$A,
                          W = datasets$df_aw[, -c(1:2)], nfolds = nfolds,
                          family = family, Qform = nice_qform,
                          gform = nice_gform)
  ests_w <- getEstimates(Y = datasets$df_w$Y, A = datasets$df_w$A,
                         W = datasets$df_w[, -c(1:2)], nfolds = nfolds,
                         family = family, Qform = nice_qform,
                         gform = nice_gform)
  return(list(est_d1 = ests_aw$ATE, est_d2 = ests_w$ATE,
              ps_sum_d1 = ests_aw$psSum, ps_sum_d2 = ests_w$psSum,
              est_ey0_d1 = ests_aw$EY0, est_ey0_d2 = ests_w$EY0,
              est_ey1_d1 = ests_aw$EY1, est_ey1_d2 = ests_w$EY1,
              est_rr_d1 = ests_aw$RR, est_rr_d2 = ests_w$RR,
              est_cate_d1 = ests_aw$CATE, est_cate_d2 = ests_w$CATE,
              cor_aw = cor_aw, cor_w = cor_w))
}


# Step 2. Draw a sample 
# @param B the number of monte-carlo replications
# @param n the sample size
# @param nfolds the number of folds for cross-validated methods
# @param seed the random number seed
# @param original_data the original dataset
# @param ps_model the treatment propensity model fit to the original dataset
# @param outcome_model the outcome regression model fit to the original dataset
runSim <- function(B=1000, n = 100, nfolds = 20, 
                   seed = 100, original_data = NULL, ps_model = NULL, 
                   outcome_model = NULL){
  # get the true treatment effect using the original fitted outcome model
  # note that this requires Y,A to be the first two columns of the dataset
  ey1 <- mean(pred.Y(A = rep(1, nrow(original_data)), 
                     W = original_data[,-(1:2)],
                     outcome_model = outcome_model))
  ey0 <- mean(pred.Y(A = rep(0, nrow(original_data)), 
                     W = original_data[,-(1:2)],
                     outcome_model = outcome_model))
  psi0.sim <- ey1 - ey0
  RR0 <- ey1 / ey0
  CATE0 <- coefficients(outcome_model)[2]
      
  
  # set up matrices to hold results
  est.d1.ATE <- matrix(nrow = B, ncol = 11)
  colnames(est.d1.ATE) <- c("unadj", "Match", "IPTW", "IPTW.unbd", "TMLE", 
                            "TMLE.unbd", "glmCor", "glmPS", "PSstrat5",  
                            "PSstrat10", "PSstrat20")
  est.d2.ATE <-  est.d1.ATE
  psSum.d1 <- matrix(nrow = B, ncol = 6)
  colnames(psSum.d1) <- c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.")
  psSum.d2 <- psSum.d1
  
  # set up parallelization and random number stream
  future::plan(multisession)
  options(future.globals.maxSize = +Inf)
  set.seed(seed)
  seeds <- future_lapply(as.list(seq_len(B)), FUN = function(x) .Random.seed,
                         future.chunk.size = Inf, future.seed = seed)
  output <- future_lapply(
    X = as.list(1:B), FUN = function(b) {
      tryCatch(
        one_sim(b = b, original_covariates = original_data[, -c(1:2)],
                n = n, original_tx = original_data$A, ps_model = ps_model,
                outcome_model = outcome_model, family = "binomial"),
        error = function(e) {
          return(list(est_d1 = rep(NA, 11), est_d2 = rep(NA, 11), 
                      ps_sum_d1 = rep(NA, 6), ps_sum_d2 = rep(NA, 6)))
        }
      )
    }, future.seed = seeds, future.stdout = structure(TRUE, drop = TRUE),
    future.conditions = "message"
  )
  est.d1.ATE <- do.call(rbind, lapply(output, function(x) x$est_d1))
  est.d2.ATE <- do.call(rbind, lapply(output, function(x) x$est_d2))
  est.d1.RR <- do.call(rbind, lapply(output, function(x) x$est_rr_d1))
  est.d2.RR <- do.call(rbind, lapply(output, function(x) x$est_rr_d2))
  est.d1.EY0 <- do.call(rbind, lapply(output, function(x) x$est_ey0_d1))
  est.d2.EY0 <- do.call(rbind, lapply(output, function(x) x$est_ey0_d2))
  est.d1.EY1 <- do.call(rbind, lapply(output, function(x) x$est_ey1_d1))
  est.d2.EY1 <- do.call(rbind, lapply(output, function(x) x$est_ey1_d2))
  est.d1.CATE <- do.call(rbind, lapply(output, function(x) x$est_cate_d1))
  est.d2.CATE <- do.call(rbind, lapply(output, function(x) x$est_cate_d2))
  
  psSum.d1 <- do.call(rbind, lapply(output, function(x) x$ps_sum_d1))
  psSum.d2 <- do.call(rbind, lapply(output, function(x) x$ps_sum_d2))
  
  cor_aw <- do.call(rbind, lapply(output, function(x) x$cor_aw))
  cor_w <- do.call(rbind, lapply(output, function(x) x$cor_w))
  
  # evaluate results
  # evaluate the true injected signal	
  res.d1.ATE <- calcResults(est.d1.ATE, psi0.sim)
  res.d2.ATE <- calcResults(est.d2.ATE, psi0.sim)	
  res.d1.RR <- calcResults(est.d1.RR, RR0)
  res.d2.RR <- calcResults(est.d2.RR, RR0)
  res.d1.EY0 <- calcResults(est.d1.EY0, ey0)
  res.d2.EY0 <- calcResults(est.d2.EY0, ey0)
  res.d1.EY1 <- calcResults(est.d1.EY1, ey1)
  res.d2.EY1 <- calcResults(est.d2.EY1, ey1)
  res.d1.CATE <- calcResults(est.d1.CATE, CATE0)
  res.d2.CATE <- calcResults(est.d2.CATE, CATE0)
  
  res.sim <- list(n = n, B = B, psSum.d1 = psSum.d1, psSum.d2 = psSum.d2,
                  cor_aw = cor_aw, cor_w = cor_w,
                  ATE= list(est.d1 = est.d1.ATE, est.d2 = est.d2.ATE,  
                            res.d1=res.d1.ATE, res.d2=res.d2.ATE,psi0 = psi0.sim),
                  RR = list(est.d1 = est.d1.RR, est.d2 = est.d2.RR,
                            res.d1 = res.d1.RR, res.d2 = res.d2.RR, psi0 = RR0),
                  EY0 = list(est.d1 = est.d1.EY0, est.d2 = est.d2.EY0,
                            res.d1 = res.d1.EY0, res.d2 = res.d2.EY0, psi0 = ey0),
                  EY1 = list(est.d1 = est.d1.EY1, est.d2 = est.d2.EY1,
                            res.d1 = res.d1.EY1, res.d2 = res.d2.EY1, psi0 = ey1),
                  CATE = list(est.d1 = est.d1.CATE, est.d2 = est.d2.CATE,
                            res.d1 = res.d1.CATE, res.d2 = res.d2.CATE, psi0 = CATE0)) 
  return(res.sim)				
}
