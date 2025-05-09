# compile results from plasmode paper simulations

# load required functions and packages ----------------------------------------
library("dplyr")
library("tibble")
library("tidyr")
library("data.table")
library("ggplot2")
library("cowplot")
theme_set(theme_cowplot())

# the following lines must be edited for your filesystem
main_dir <- "<replace with your filesystem root>"
proj_root <- paste0(main_dir, "<replace with path to your project root>")
setwd(proj_root)
production_code_dir <- paste0(proj_root, "R code/")
production_output_dir <- paste0(proj_root, "Results/plasmode_aw/")
data_dir <- paste0(proj_root, "Data/")

source(paste0(production_code_dir, "00_plasmode_utils.R"))

# read in results -------------------------------------------------------------
outcome_rates <- c(0.05, 0.15)
ns <- 10000
B <- 1e5
n_rate <- expand.grid(n = ns, rate = outcome_rates)
output_dirs <- lapply(1:nrow(n_rate), function(i) {
  paste0(production_output_dir, "n", n_rate$n[i], "_rate", n_rate$rate[i], "/")
})

full_results <- lapply(output_dirs, function(x) {
  readRDS(paste0(x, "results_b", B, ".rds"))
})

# est.d1 is sampling (A,W) together
all_ests_aw_ATE <- data.table::rbindlist(
  lapply(1:length(full_results), function(i) {
    data.table::rbindlist(
      lapply(1:length(full_results[[i]]), function(j) {
        data.frame(n = n_rate$n[i], rate = n_rate$rate[i], 
                   full_results[[i]][[j]]$ATE$est.d1)
      })
    )
  })
)
all_ests_aw_RR <- data.table::rbindlist(
  lapply(1:length(full_results), function(i) {
    data.table::rbindlist(
      lapply(1:length(full_results[[i]]), function(j) {
        data.frame(n = n_rate$n[i], rate = n_rate$rate[i], 
                   full_results[[i]][[j]]$RR$est.d1)
      })
    )
  })
)
all_ests_aw_EY0 <- data.table::rbindlist(
  lapply(1:length(full_results), function(i) {
    data.table::rbindlist(
      lapply(1:length(full_results[[i]]), function(j) {
        data.frame(n = n_rate$n[i], rate = n_rate$rate[i], 
                   full_results[[i]][[j]]$EY0$est.d1)
      })
    )
  })
)
all_ests_aw_EY1 <- data.table::rbindlist(
  lapply(1:length(full_results), function(i) {
    data.table::rbindlist(
      lapply(1:length(full_results[[i]]), function(j) {
        data.frame(n = n_rate$n[i], rate = n_rate$rate[i], 
                   full_results[[i]][[j]]$EY1$est.d1)
      })
    )
  })
)
all_ests_aw_CATE <- data.table::rbindlist(
  lapply(1:length(full_results), function(i) {
    data.table::rbindlist(
      lapply(1:length(full_results[[i]]), function(j) {
        data.frame(n = n_rate$n[i], rate = n_rate$rate[i], 
                   full_results[[i]][[j]]$CATE$est.d1)
      })
    )
  })
)
all_ests_aw_cor_aw <- data.table::rbindlist(
  lapply(1:length(full_results), function(i) {
    data.table::rbindlist(
      lapply(1:length(full_results[[i]]), function(j) {
        data.frame(n = n_rate$n[i], rate = n_rate$rate[i], 
                   full_results[[i]][[j]]$cor_aw)
      })
    )
  })
)


# est.d2 is sampling W, then generating A
all_ests_w_ATE <- data.table::rbindlist(
  lapply(1:length(full_results), function(i) {
    data.table::rbindlist(
      lapply(1:length(full_results[[i]]), function(j) {
        data.frame(n = n_rate$n[i], rate = n_rate$rate[i], 
                   full_results[[i]][[j]]$ATE$est.d2)
      })
    )
  })
)
all_ests_w_RR <- data.table::rbindlist(
  lapply(1:length(full_results), function(i) {
    data.table::rbindlist(
      lapply(1:length(full_results[[i]]), function(j) {
        data.frame(n = n_rate$n[i], rate = n_rate$rate[i], 
                   full_results[[i]][[j]]$RR$est.d2)
      })
    )
  })
)
all_ests_w_EY0 <- data.table::rbindlist(
  lapply(1:length(full_results), function(i) {
    data.table::rbindlist(
      lapply(1:length(full_results[[i]]), function(j) {
        data.frame(n = n_rate$n[i], rate = n_rate$rate[i], 
                   full_results[[i]][[j]]$EY0$est.d2)
      })
    )
  })
)
all_ests_w_EY1 <- data.table::rbindlist(
  lapply(1:length(full_results), function(i) {
    data.table::rbindlist(
      lapply(1:length(full_results[[i]]), function(j) {
        data.frame(n = n_rate$n[i], rate = n_rate$rate[i], 
                   full_results[[i]][[j]]$EY1$est.d2)
      })
    )
  })
)
all_ests_w_CATE <- data.table::rbindlist(
  lapply(1:length(full_results), function(i) {
    data.table::rbindlist(
      lapply(1:length(full_results[[i]]), function(j) {
        data.frame(n = n_rate$n[i], rate = n_rate$rate[i], 
                   full_results[[i]][[j]]$CATE$est.d2)
      })
    )
  })
)
all_ests_w_cor_w <- data.table::rbindlist(
  lapply(1:length(full_results), function(i) {
    data.table::rbindlist(
      lapply(1:length(full_results[[i]]), function(j) {
        data.frame(n = n_rate$n[i], rate = n_rate$rate[i], 
                   full_results[[i]][[j]]$cor_w)
      })
    )
  })
)

# truths
psi0 <- tibble::tibble(
  rate = rep(outcome_rates, 5),
  estimand = c(rep("ATE", 2), rep("RR", 2),
               rep("EY1", 2), rep("EY0", 2),
               rep("CATE", 2)),
  psi0 = c(full_results[[1]][[1]]$ATE$psi0, full_results[[2]][[1]]$ATE$psi0,
           full_results[[1]][[1]]$RR$psi0, full_results[[2]][[1]]$RR$psi0,
           full_results[[1]][[1]]$EY1$psi0, full_results[[2]][[1]]$EY1$psi0,
           full_results[[1]][[1]]$EY0$psi0, full_results[[2]][[1]]$EY0$psi0,
           full_results[[1]][[1]]$CATE$psi0, full_results[[2]][[1]]$CATE$psi0)
)

all_ests <- dplyr::bind_rows(
  all_ests_aw_ATE %>% mutate(sampling = "AW", estimand = "ATE") %>% 
    rename_with(~stringr::str_remove(., '.A')) %>% 
    select(-contains("PSstrat")),
  all_ests_w_ATE %>% mutate(sampling = "W", estimand = "ATE") %>% 
    rename_with(~stringr::str_remove(., '.A')) %>% 
    select(-contains("PSstrat")),
  all_ests_aw_RR %>% mutate(sampling = "AW", estimand = "RR"),
  all_ests_w_RR %>% mutate(sampling = "W", estimand = "RR"),
  all_ests_aw_EY1 %>% mutate(sampling = "AW", estimand = "EY1"),
  all_ests_w_EY1 %>% mutate(sampling = "W", estimand = "EY1"),
  all_ests_aw_EY0 %>% mutate(sampling = "AW", estimand = "EY0"),
  all_ests_w_EY0 %>% mutate(sampling = "W", estimand = "EY0"),
  all_ests_aw_CATE %>% mutate(sampling = "AW", estimand = "CATE"),
  all_ests_w_CATE %>% mutate(sampling = "W", estimand = "CATE")
) %>% 
  left_join(psi0, by = c("rate", "estimand"))

all_cors <- dplyr::bind_rows(
  all_ests_aw_cor_aw %>% mutate(sampling = "Sample Treatment"),
  all_ests_w_cor_w %>% mutate(sampling = "Generate Treatment")
) %>% 
  pivot_longer(cols = sex:AgeIndexsq, 
               names_to = "Variable",
               values_to = "Correlation") %>% 
  mutate(Variable = case_when(
    Variable == "sex" ~ "Sex",
    Variable == "AgeAtIndex" ~ "Age",
    Variable == "charlson_1" ~ "Charlson = 1",
    Variable == "charlson_2" ~ "Charlson = 2",
    Variable == "charlson_3plus" ~ "Charlson = 3+",
    Variable == "aud_dx_priorYr" ~ "AUD Dx",
    Variable == "anx_dx_priorYr" ~ "ANX Dx",
    Variable == "priorSH" ~ "Prior SH",
    Variable == "MHIP_prior5yr" ~ "Prior MH Hosp",
    Variable == "phq8_6to10" ~ "PHQ8 = 6-10",
    Variable == "phq8_11to15" ~ "PHQ8 = 11-15",
    Variable == "phq8_16to20" ~ "PHQ8 = 16-20",
    Variable == "phq8_20plus" ~ "PHQ8 = 20+",
    Variable == "phq9_some" ~ "PHQ9 = 1",
    Variable == "phq9_most" ~ "PHQ9 = 2",
    Variable == "phq9_everyday" ~ "PHQ9 = 3",
    Variable == "AgeIndexsq" ~ "Age sq."
  ))

# write out tables ------------------------------------------------------------
# get summaries over the B monte-carlo replications
all_ests_long <- all_ests %>% 
  select(n, rate, sampling, estimand, psi0, unadj:est.glmPS, glm.A) %>% 
  pivot_longer(unadj:glm.A,
               names_to = "estimator",
               values_to = "estimate") %>% 
  # remove crazy CATE estimates
  filter(!(grepl("CATE", estimand) & abs(estimate) > 10)) %>% 
  group_by(n, rate, sampling, estimand, estimator, psi0) %>% 
  filter(!(is.na(estimate) & (estimator == "glm.A") & (estimand != "CATE"))) %>% 
  filter(!(is.na(estimate) & (estimator != "glm.A") & (estimand == "CATE"))) %>% 
  mutate(estimand = ifelse(estimand == "CATE", "logcOR", estimand))

ses <- all_ests_long %>% 
  group_by(n, rate, sampling, estimand, estimator, psi0) %>% 
  summarize(SE = sd(estimate), .groups = "drop")

summs <- all_ests_long %>% 
  left_join(ses, by = c("n", "rate", "sampling", "estimand", "estimator", "psi0")) %>% 
  group_by(n, rate, sampling, estimand, estimator, psi0, SE) %>% 
  summarize(`Mean Bias` = mean(estimate - psi0),
            `Median Bias` = median(estimate - psi0),
            MSE = mean((estimate - psi0) ^ 2), 
            `Cov.` = mean(estimate - qnorm(0.975) * SE <= psi0 & psi0 <= estimate + qnorm(0.975) * SE), 
            .groups = "drop") %>% 
  mutate(`% Bias` = `Mean Bias` / psi0 * 100, 
         `bias:SE` = abs(`Mean Bias`) / SE,
         RMSE = sqrt(MSE),
         `Cov.` = `Cov.` * 100,
         estimator = factor(estimator, 
                            levels = c("unadj", "match", "iptw",
                                       "iptw.unbd", "tmle", "tmle.unbd",
                                       "est.glm", "est.glmPS", "glm.A"),
                            labels = c("Unadj", "Match",
                                       "IPTW", "IPTW (unbounded)",
                                       "TMLE", "TMLE (unbounded)",
                                       "glmCM", "glmPS",
                                       "glmCM")),
         num_estimator = case_when(
           estimator == "Unadj" ~ 1,
           estimator == "Match" ~ 2,
           estimator == "IPTW" ~ 3,
           estimator == "IPTW (unbounded)" ~ 4,
           estimator == "TMLE" ~ 5,
           estimator == "TMLE (unbounded)" ~ 6,
           estimator == "glmCM" ~ 7,
           estimator == "glmPS" ~ 8
         ),
         num_estimand = case_when(
           estimand == "ATE" ~ 1,
           estimand == "RR" ~ 2,
           estimand == "EY0" ~ 3,
           estimand == "EY1" ~ 4,
           estimand == "logcOR" ~ 5
         ))

# full results table -----------------------------------------------------------
# for each outcome and estimand, generate a full results table
# edit 20250127: put both ATE and RR in the same table
abbreviation_txt_prefix <- paste0("Abbreviations: ",
                           "KPWA: Kaiser Permanente Washington; ")
ate_rr_txt <- "ATE: average treatment effect; RR: relative risk; "
other_estimand_txt <- paste0("EY0: expected value of outcome under no treatment; ",
                             "EY1: expected value of outcome under treatment",
                             "logcOR: conditional log odds ratio; ")
abbreviation_txt_suffix <- paste0(
                           "SE: standard error; `RMSE': root mean squared error; ",
                           "`Bias:SE': ratio of bias to standard error; ",
                           "`Cov.': confidence interval coverage; ",
                           "Unadj: unadjusted; ",
                           "Match: propensity score matching; ",
                           "IPTW: inverse probability of treatment weighting; ",
                           "TMLE: targeted maximum likelihood estimation; ",
                           "glmCM: generalized linear model, correctly specified; ",
                           "glmPS: generalized linear model, adjusted for propensity score."
                           )
one_over_rootn <- 1 / sqrt(summs$n[1])
n_mc_txt <- paste0(" Results for n = 10,000 ($1/\\\\sqrt{n}=$", one_over_rootn, ") with 100,000 Monte-Carlo iterations. ")
other_estimands <- c("EY0", "EY1", "logcOR")
font_size <- 11

for (i in 1:length(unique(summs$rate))) {
  this_rate <- unique(summs$rate)[i]
  # ATE
  this_summ_ate <- summs %>% 
    filter(rate == this_rate, estimand == "ATE") %>% 
    select(-n, -rate, -estimand) 
  this_psi0_ate <- round(this_summ_ate$psi0[1], 3)
  wide_summ_ate <- this_summ_ate %>% 
    select(-psi0) %>% 
    pivot_wider(id_cols = c(estimator, num_estimator),
                names_from = sampling,
                values_from = c(`Mean Bias`, `Median Bias`, SE, MSE, RMSE, `% Bias`, `bias:SE`, `Cov.`),
                names_glue = "{.value} ({sampling})") %>% 
    arrange(num_estimator) %>% 
    select(-num_estimator) %>% 
    select(estimator, contains("AW"), !contains("AW")) %>% 
    rename(Estimator = estimator)
  sample_aw_ate <- wide_summ_ate %>% 
    select(ends_with("(AW)")) %>% 
    rename_with(~ gsub(" (AW)", "", .x, fixed = TRUE)) %>% 
    mutate(across(`Mean Bias`:`bias:SE`, ~ ifelse(abs(.x) < .001, "$<$ 0.001", 
                                                  sprintf("%.3f", .x))))
  sample_w_ate <- wide_summ_ate %>% 
    select(ends_with("(W)")) %>% 
    rename_with(~ gsub(" (W)", "", .x, fixed = TRUE)) %>% 
    mutate(across(`Mean Bias`:`bias:SE`, ~ ifelse(abs(.x) < .001, "$<$ 0.001", sprintf("%.3f", .x))))
  wide_summ2_ate <- data.frame(
    Estimand = "ATE",
    Estimator = wide_summ_ate %>% 
      pull(Estimator),
    sample_aw_ate,
    sample_w_ate
  )
  # RR
  this_summ_rr <- summs %>% 
    filter(rate == this_rate, estimand == "RR") %>% 
    select(-n, -rate, -estimand) 
  this_psi0_rr <- round(this_summ_rr$psi0[1], 3)
  wide_summ_rr <- this_summ_rr %>% 
    select(-psi0) %>% 
    pivot_wider(id_cols = c(estimator, num_estimator),
                names_from = sampling,
                values_from = c(`Mean Bias`, `Median Bias`, SE, MSE, RMSE, `% Bias`, `bias:SE`, `Cov.`),
                names_glue = "{.value} ({sampling})") %>% 
    arrange(num_estimator) %>% 
    select(-num_estimator) %>% 
    select(estimator, contains("AW"), !contains("AW")) %>% 
    rename(Estimator = estimator)
  sample_aw_rr <- wide_summ_rr %>% 
    select(ends_with("(AW)")) %>% 
    rename_with(~ gsub(" (AW)", "", .x, fixed = TRUE)) %>% 
    mutate(across(`Mean Bias`:`bias:SE`, ~ ifelse(abs(.x) < .001, "$<$ 0.001", sprintf("%.3f", .x))))
  sample_w_rr <- wide_summ_rr %>% 
    select(ends_with("(W)")) %>% 
    rename_with(~ gsub(" (W)", "", .x, fixed = TRUE)) %>% 
    mutate(across(`Mean Bias`:`bias:SE`, ~ ifelse(abs(.x) < .001, "$<$ 0.001", sprintf("%.3f", .x))))
  wide_summ2_rr <- data.frame(
    Estimand = "RR",
    Estimator = wide_summ_rr %>% 
      pull(Estimator),
    sample_aw_rr,
    sample_w_rr
  )
  wide_summ2 <- rbind(wide_summ2_ate, wide_summ2_rr)
  
  nice_names <- names(sample_aw_ate)
  nice_names[6] <- "\\% Bias"
  
  true_value_txt <- paste0("The true value of the ATE is ", this_psi0_ate, "; the true value of the RR is ", this_psi0_rr, ".")
  
  # reduced table, for main paper
  wide_summ2 %>% 
    filter(!grepl("unbounded", Estimator)) %>% 
    select(Estimand, Estimator, X..Bias, SE, RMSE, bias.SE, `Cov.`, X..Bias.1, SE.1, RMSE.1, bias.SE.1, Cov..1) %>% 
    knitr::kable(format = "latex", digits = 1, escape = FALSE,
                 booktabs = TRUE, linesep = "", caption = paste0(
                   "Results of the KPWA-based plasmode simulations for the ",
                   this_rate * 100, "\\% outcome (ATE and RR).\\label{tab:kpwa_plasmode_r", this_rate, "}"
                 ),
                 col.names = c("Estimand", "Estimator",
                               rep(c("\\% Bias", "SE", "RMSE", "bias:SE", "Cov."), 2)
                 )) %>% 
    kableExtra::kable_styling(font_size = font_size, latex_options = "scale_down") %>% 
    kableExtra::add_header_above(c(" " = 2, "Sample Treatment" = 5, "Generate Treatment" = 5)) %>%
    kableExtra::footnote(general = paste0(abbreviation_txt_prefix, ate_rr_txt, abbreviation_txt_suffix, n_mc_txt, true_value_txt),
                         general_title = "",
                         threeparttable = TRUE, escape = FALSE) %>% 
    kableExtra::collapse_rows(columns = 1) %>% 
    kableExtra::save_kable(file = paste0(production_output_dir, "kpwa_plasmode_results_r", this_rate, ".tex"))
  # other estimands for supplement
  this_summ_other <- summs %>%
    filter(rate == this_rate, estimand %in% other_estimands) %>% 
    select(-n, -rate)
  these_psi0s_other <- this_summ_other %>% 
    group_by(estimand) %>% 
    slice(1) %>% 
    select(psi0) %>% 
    mutate(psi0 = round(psi0, 3))
  psi0_other_txt <- paste0("The true value of EY0 is ", these_psi0s_other$psi0[these_psi0s_other$estimand == "EY0"],
                           "; the true value of EY1 is ", these_psi0s_other$psi0[these_psi0s_other$estimand == "EY1"],
                           "; the true value of logcOR is ", these_psi0s_other$psi0[these_psi0s_other$estimand == "logcOR"], ".")
  wide_summ_other <- this_summ_other %>% 
    select(-psi0) %>% 
    pivot_wider(id_cols = c(estimand, num_estimand, estimator, num_estimator),
                names_from = sampling,
                values_from = c(`Mean Bias`, `Median Bias`, SE, MSE, RMSE, `% Bias`, `bias:SE`, `Cov.`),
                names_glue = "{.value} ({sampling})") %>% 
    arrange(num_estimand, num_estimator) %>% 
    select(-num_estimand, -num_estimator) %>% 
    select(estimand, estimator, contains("AW"), !contains("AW")) %>% 
    rename(Estimand = estimand, Estimator = estimator)
  sample_aw_other <- wide_summ_other %>% 
    select(ends_with("(AW)")) %>% 
    rename_with(~ gsub(" (AW)", "", .x, fixed = TRUE)) %>% 
    mutate(across(`Mean Bias`:`bias:SE`, ~ ifelse(abs(.x) < .001, "$<$ 0.001", sprintf("%.3f", .x))))
  sample_w_other <- wide_summ_other %>% 
    select(ends_with("(W)")) %>% 
    rename_with(~ gsub(" (W)", "", .x, fixed = TRUE)) %>% 
    mutate(across(`Mean Bias`:`bias:SE`, ~ ifelse(abs(.x) < .001, "$<$ 0.001", sprintf("%.3f", .x))))
  wide_summ2_other <- data.frame(
    Estimand = wide_summ_other %>% 
      pull(Estimand),
    Estimator = wide_summ_other %>% 
      pull(Estimator),
    sample_aw_other,
    sample_w_other
  )
  wide_summ2_other %>% 
    filter(!grepl("unbounded", Estimator)) %>% 
    select(Estimand, Estimator, X..Bias, SE, RMSE, bias.SE, Cov., X..Bias.1, SE.1, RMSE.1, bias.SE.1, Cov..1) %>% 
    knitr::kable(format = "latex", digits = 1, escape = FALSE,
                 booktabs = TRUE, linesep = "", caption = paste0(
                   "Results of the KPWA-based plasmode simulations for the ",
                   this_rate * 100, "\\% outcome.\\label{tab:kpwa_plasmode_r", this_rate, "_supp}"
                 ),
                 col.names = c("Estimand", "Estimator",
                               rep(c("\\% Bias", "SE", "RMSE", "bias:SE", "Coverage"), 2)
                 )) %>% 
    kableExtra::kable_styling(font_size = font_size, latex_options = "scale_down") %>% 
    kableExtra::add_header_above(c(" " = 2, "Sample Treatment" = 5, "Generate Treatment" = 5)) %>%
    kableExtra::footnote(general = paste0(abbreviation_txt_prefix, other_estimand_txt, 
                                          abbreviation_txt_suffix, n_mc_txt, psi0_other_txt),
                         general_title = "",
                         threeparttable = TRUE, escape = FALSE) %>% 
    kableExtra::collapse_rows(columns = 1) %>% 
    kableExtra::save_kable(file = paste0(production_output_dir, "kpwa_plasmode_results_r", this_rate, "_supp.tex"))
}

# bias to SE with all outcomes in a table --------------------------------------
for (i in 1:length(unique(summs$estimand))) {
  this_estimand <- unique(summs$estimand)[i]
  this_summ <- summs %>% 
    filter(estimand == this_estimand)
  nice_bias_to_se <- this_summ %>% 
    select(n, rate, estimator, sampling, `bias:SE`, num_estimator) %>% 
    mutate(sampling = ifelse(grepl("AW", sampling), "Sample Treatment", "Generate Treatment")) %>% 
    pivot_wider(id_cols = c(n, rate, estimator, num_estimator),
                names_from = sampling,
                values_from = `bias:SE`,
                names_glue = "bias:SE ({sampling})") %>% 
    arrange(rate, num_estimator) %>% 
    select(-num_estimator)
  write.csv(
    nice_bias_to_se,
    paste0(production_output_dir, "bias_to_se_", this_estimand, ".csv")
  ) 
}

# plot of correlations ---------------------------------------------------------
cor_plot <- all_cors %>% 
  ggplot(aes(y = Correlation, x = Variable, color = sampling)) +
  geom_boxplot() +
  scale_color_manual(values = c("gray40", "gray70")) +
  labs(color = "Approach") + 
  ylab("Correlation with treatment assignment") +
  facet_grid(rows = vars(rate)) +
  scale_x_discrete(guide = guide_axis(angle = 60))
ggsave(file = paste0(production_output_dir, "kpwa_plasmode_correlations.png"),
       cor_plot, width = 10, height = 10, units = "in")
