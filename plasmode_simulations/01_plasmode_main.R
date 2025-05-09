# plasmode simulation comparing sampling A (treatment) and W (confounders) together
# versus sampling W and generating A

# load required packages and functions -----------------------------------------
library("future.apply")
library("optparse")
library("MatchIt")
library("tmle")
library("dplyr")

# the following lines must be edited for your filesystem
main_dir <- "<replace with your filesystem root>"
proj_root <- paste0(main_dir, "<replace with path to your project root>")
setwd(proj_root)
production_code_dir <- paste0(proj_root, "R code/")
production_output_dir <- paste0(proj_root, "Results/")
data_dir <- paste0(proj_root, "Data/")

source(paste0(production_code_dir, "00_plasmode_utils.R"))

parser <- OptionParser()
parser <- add_option(parser, "--nreps-total", default = 2500, help = "Number of simulation replicates in total")
parser <- add_option(parser, "--nreps-per-job", default = 2500 / 10, help = "Number of simulation replicates to run at once before saving")
parser <- add_option(parser, "--n", default = 1e4, help = "Sample size")
parser <- add_option(parser, "--outcome-rate", default = 0.15, help = "Outcome rate")
args <- parse_args(parser, convert_hyphens_to_underscores = TRUE)

num_jobs <- args$nreps_total / args$nreps_per_job

this_output_dir <- paste0(production_output_dir, "plasmode_aw/", "n", args$n, "_rate", args$outcome_rate, "/")
if (!dir.exists(this_output_dir)) {
  dir.create(this_output_dir, recursive = TRUE)
}

# get regression parameters for propensity score, outcome models ---------------------------
# this loads in the original dataset, called 'data'
load(file = paste0(data_dir, "plasmode data sets/", "AD_PT_plasmode_full.Rdata"))

# rename numEpiType to be A, SH_HOSP_1826days to be Y; convert categorical/factor variables to dummy variables
data <- data %>% 
  rename(A = numEpiType, Y = SH_HOSP_1826days) %>% 
  mutate(charlson_1 = as.numeric(Charlson_cat == "1"),
         charlson_2 = as.numeric(Charlson_cat == "2"),
         charlson_3plus = as.numeric(Charlson_cat == "3+"),
         phq8_6to10 = as.numeric(IndexTrtPHQ8_score_cat == "6-10"),
         phq8_11to15 = as.numeric(IndexTrtPHQ8_score_cat == "11-15"),
         phq8_16to20 = as.numeric(IndexTrtPHQ8_score_cat == "16-20"),
         phq8_20plus = as.numeric(IndexTrtPHQ8_score_cat == "20+"),
         phq9_some = as.numeric(IndexTrtPHQ_item9_score_cat == "Some of the days"),
         phq9_most = as.numeric(IndexTrtPHQ_item9_score_cat == "Most of the days"),
         phq9_everyday = as.numeric(IndexTrtPHQ_item9_score_cat == "Nearly everyday")) %>% 
  select(-contains("_cat"))

# propensity of treatment
m.A <- glm(A ~ sex + AgeAtIndex + AgeIndexsq + 
             charlson_1 + charlson_2 + charlson_3plus +
             anx_dx_priorYr + (charlson_1 + charlson_2 + charlson_3plus) * anx_dx_priorYr + 
             aud_dx_priorYr + priorSH + MHIP_prior5yr + sex * AgeAtIndex +
             priorSH * sex + priorSH * AgeAtIndex +
             (charlson_1 + charlson_2 + charlson_3plus) * AgeAtIndex + 
             phq8_6to10 + phq8_11to15 + phq8_16to20 + phq8_20plus +
             phq9_some + phq9_most + phq9_everyday +
             (phq9_some + phq9_most + phq9_everyday) * sex +
             (phq9_some + phq9_most + phq9_everyday) * priorSH, 
           data = data, family = "binomial")
# outcome regression  
m.Y <- glm(Y ~ A + sex + AgeAtIndex +
             charlson_1 + charlson_2 + charlson_3plus + aud_dx_priorYr + 
             anx_dx_priorYr + priorSH + MHIP_prior5yr + 
             phq8_6to10 + phq8_11to15 + phq8_16to20 + phq8_20plus +
             phq9_some + phq9_most + phq9_everyday + AgeIndexsq +
             (charlson_1 + charlson_2 + charlson_3plus) * anx_dx_priorYr + MHIP_prior5yr +
             sex * AgeAtIndex + priorSH * sex + priorSH * AgeAtIndex +
             (charlson_1 + charlson_2 + charlson_3plus) * AgeAtIndex + 
             (phq9_some + phq9_most + phq9_everyday) * sex +
             (phq9_some + phq9_most + phq9_everyday) * priorSH,
           data = data, family = "binomial")

W <- data %>%
  select(sex, AgeAtIndex, starts_with("charlson_"), aud_dx_priorYr, anx_dx_priorYr,
         priorSH, MHIP_prior5yr, 
         starts_with("phq8_"), starts_with("phq9_"), AgeIndexsq)
ey0 <- mean(predict(m.Y, newdata = data.frame(A = 0, W), type = "response"))
ey1 <- mean(predict(m.Y, newdata = data.frame(A = 1, W), type = "response"))
psi0_orig <- ey1 - ey0
# new coefficients for outcome model with rate 15%
m.Y2 <- m.Y
m.Y2$coefficients[c(1, 2)] <- c(-1.32, -1)
m.Y3 <- m.Y2
# outcome rate 5%
m.Y3$coefficients[c(1,2)] <- c(-2.35, -3.1)

if (args$outcome_rate == 0.15) {
  # generate a new "original" dataset
  set.seed(20241121)
  original_data_list <- gen_plasmode_dataset(original_covariates = W, n = args$n, 
                                             original_tx = data$A, 
                                             ps_model = m.A, outcome_model = m.Y2)
} else {
  set.seed(20241126)
  original_data_list <- gen_plasmode_dataset(original_covariates = W, n = args$n, 
                                             original_tx = data$A, 
                                             ps_model = m.A, outcome_model = m.Y3)
}
original_data <- original_data_list$df_w
# fit a new propensity score model and outcome regression model
ps_model <- glm(A ~ sex + AgeAtIndex + AgeIndexsq + 
                  charlson_1 + charlson_2 + charlson_3plus +
                  anx_dx_priorYr + (charlson_1 + charlson_2 + charlson_3plus) * anx_dx_priorYr + 
                  aud_dx_priorYr + priorSH + MHIP_prior5yr + sex * AgeAtIndex +
                  priorSH * sex + priorSH * AgeAtIndex +
                  (charlson_1 + charlson_2 + charlson_3plus) * AgeAtIndex + 
                  phq8_6to10 + phq8_11to15 + phq8_16to20 + phq8_20plus +
                  phq9_some + phq9_most + phq9_everyday +
                  (phq9_some + phq9_most + phq9_everyday) * sex +
                  (phq9_some + phq9_most + phq9_everyday) * priorSH, 
                data = original_data, family = "binomial")

outcome_model <- glm(Y ~ A + sex + AgeAtIndex +
                       charlson_1 + charlson_2 + charlson_3plus + aud_dx_priorYr + 
                       anx_dx_priorYr + priorSH + MHIP_prior5yr + 
                       phq8_6to10 + phq8_11to15 + phq8_16to20 + phq8_20plus +
                       phq9_some + phq9_most + phq9_everyday + AgeIndexsq +
                       (charlson_1 + charlson_2 + charlson_3plus) * anx_dx_priorYr + MHIP_prior5yr +
                       sex * AgeAtIndex + priorSH * sex + priorSH * AgeAtIndex +
                       (charlson_1 + charlson_2 + charlson_3plus) * AgeAtIndex + 
                       (phq9_some + phq9_most + phq9_everyday) * sex +
                       (phq9_some + phq9_most + phq9_everyday) * priorSH,,
                     data = original_data, family = "binomial")
ey0_plasmode <- mean(predict(outcome_model, newdata = data.frame(A = 0, original_data[, -c(1:2)]), type = "response"))
ey1_plasmode <- mean(predict(outcome_model, newdata = data.frame(A = 1, original_data[, -c(1:2)]), type = "response"))
psi0_plasmode <- ey1_plasmode - ey0_plasmode
cat("Original data ATE:\n")
psi0_orig
cat("New plasmode data ATE,", args$outcome_rate * 100, "%:\n")
psi0_plasmode
cat("Original data RR:\n")
ey1 / ey0
cat("New plasmode RR:\n")
ey1_plasmode / ey0_plasmode
# run the simulation -----------------------------------------------------------
# run the simulation, saving output every 1000 replications
cat("\nRunning simulation\n")
results <- vector("list", length = num_jobs)
set.seed(20241121)
seeds <- round(runif(num_jobs, 1e4, 1e5))
start <- Sys.time()
for (i in 1:num_jobs) {
  results[[i]] <- runSim(B = args$nreps_per_job, n = args$n,
                         original_data = original_data,
                         ps_model = ps_model, outcome_model = outcome_model,  
                         nfolds = 20, seed = seeds[i])
  saveRDS(results, file = paste0(this_output_dir, "results_b", args$nreps_total, ".rds"))
}
end <- Sys.time()
cat("\nElapsed time: ", format(end - start), "\n")

# save results -----------------------------------------------------------------
saveRDS(results, file = paste0(this_output_dir, "results_b", args$nreps_total, ".rds"))