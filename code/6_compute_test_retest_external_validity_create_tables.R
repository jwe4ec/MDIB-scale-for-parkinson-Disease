# ---------------------------------------------------------------------------- #
# Compute Test-Retest Reliability and External Validity and Create Tables
# Original code by: Jeremy W. Eberle
# Minor modifications by: Dandan Tang
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
# Notes ----
# ---------------------------------------------------------------------------- #

# Before running script, restart R (CTRL+SHIFT+F10 on Windows) and set working 
# directory to parent folder

# ---------------------------------------------------------------------------- #
# Check correct R version and load packages ----
# ---------------------------------------------------------------------------- #

# Load custom functions

source("~/code/1a_define_functions.R")

# Check correct R version, load groundhog package, and specify groundhog_day

groundhog_day <- version_control()

# Load packages

pkgs <- c("cocor", "rockchalk", "r2glmm", "flextable", "officer", "ftExtra")
groundhog.library(pkgs, groundhog_day)

# Set "flextable" package defaults and load "officer" package properties

source("~/code/1b_set_flextable_defaults.R")
source("~/code/1c_set_officer_properties.R")

# Set seed

set.seed(1234)

# ---------------------------------------------------------------------------- #
# Define functions used in script ----
# ---------------------------------------------------------------------------- #

# Define function to compute Pearson product-moment correlation between one
# variable ("x_var") and one or more others ("y_vars") from wide data frame

cor_test_one_to_many <- function(dat, x_var, y_vars, metric) {
  r     <- vector("double", length(y_vars))
  t     <- vector("double", length(y_vars))
  df    <- vector("double", length(y_vars))
  p     <- vector("double", length(y_vars))
  ci_ll <- vector("double", length(y_vars))
  ci_ul <- vector("double", length(y_vars))
  n     <- vector("double", length(y_vars))
  
  for (i in 1:length(y_vars)) {
    cor_test <- cor.test(dat[, x_var],
                         dat[, y_vars[i]],
                         type = "pearson")
    
    r[i]     <- round(cor_test$estimate, 2)
    t[i]     <- round(cor_test$statistic, 2)
    df[i]    <- cor_test$parameter
    p[i]     <- round(cor_test$p.value, 3)
    ci_ll[i] <- round(cor_test$conf.int[1], 2)
    ci_ul[i] <- round(cor_test$conf.int[2], 2)
    
    n[i]     <- sum(complete.cases(dat[c(x_var, y_vars[i])]))
  }
  
  res <- data.frame(metric, x_var, y_var = y_vars, 
                    r, t, df, p, ci_ll, ci_ul, n)
  
  # Round selected columns and ensure two decimal digits are printed
  
  two_digit_vars   <- c("r", "t", "ci_ll", "ci_ul")
  three_digit_vars <- "p"
  
  res[, two_digit_vars]   <- format(round(res[, two_digit_vars],   2),
                                    nsmall = 2, trim = TRUE)
  res[, three_digit_vars] <- format(round(res[, three_digit_vars], 3),
                                    nsmall = 3, trim = TRUE)
  
  # Format CIs
  
  res$ci <- paste0("[", res$ci_ll, ", ", res$ci_ul, "]")
  
  # Remove redundant columns
  
  rm_cols <- c("ci_ll", "ci_ul")
  
  res <- res[, !(names(res) %in% rm_cols)]
  
  # Rearrange columns
  
  res <- res[, c("metric", "x_var", "y_var", "r", "t", "df", "p", "ci", "n")]
  
  return(res)
}

# TODO (remove function if not ultimately used): Define function to compute 
# semipartial R-squared estimates and convert to semipartial r

compute_sp_r <- function(model) {
  sp_rsq <- r2beta(model, partial = TRUE, method = "lm")
  sp_rsq <- sp_rsq[sp_rsq$Effect != "Model", ]
  
  sp_r_df <- data.frame(effect          = sp_rsq$Effect,
                        sp_rsq          = sp_rsq$Rsq,
                        sp_rsq_lower.CL = sp_rsq$lower.CL,
                        sp_rsq_upper.CL = sp_rsq$upper.CL,
                        sp_r            = sqrt(sp_rsq$Rsq),
                        sp_r_lower.CL   = sqrt(sp_rsq$lower.CL),
                        sp_r_upper.CL   = sqrt(sp_rsq$upper.CL))
  
  # Maintain predictor order from model (don't reorder in descending R-squared)
  
  sp_r_df <- sp_r_df[match(sp_r_df$effect, attr(model$terms, "term.labels")), ]
  
  sp_r_df <- rbind(c("(Intercept)", rep(NA, length(sp_r_df) - 1)),
                   sp_r_df)
  
  return(sp_r_df)
}





# Define function to run hierarchical multiple regression

run_hmr <- function(dat, outcome, bbsiq_scale, mdib_scale) {
  fml_step0 <- as.formula(paste(outcome, "~", 1))
  fml_step1 <- as.formula(paste(outcome, "~", bbsiq_scale))
  fml_step2 <- as.formula(paste(outcome, "~", bbsiq_scale, "+", mdib_scale))
  
  # Need to use "do.call" so that "fml" is evaluated before being sent to "lm"
  # but quote "dat" so it is not. Otherwise, "meanCenter" lacks the evaluated
  # "fml" call to "lm" and returns an error.
  
  lm_step0 <- do.call("lm", list(fml_step0, data = quote(dat)))
  lm_step1 <- do.call("lm", list(fml_step1, data = quote(dat)))
  lm_step2 <- do.call("lm", list(fml_step2, data = quote(dat)))
  
  lm_step1_mc <- meanCenter(lm_step1, centerOnlyInteractors = FALSE)
  lm_step2_mc <- meanCenter(lm_step2, centerOnlyInteractors = FALSE)
  
  anova_diff <- anova(lm_step1_mc, lm_step2_mc)
  
  summary_lm_step1_mc <- summary(lm_step1_mc)
  summary_lm_step2_mc <- summary(lm_step2_mc)
  
  confint_lm_step1_mc <- as.data.frame(confint(lm_step1_mc))
  confint_lm_step2_mc <- as.data.frame(confint(lm_step2_mc))
  
  # Compute zero-order r between BBSIQ score and outcome (at Step 1)
  
  bbsiq_scale_r_lm_step1 <- sqrt(summary_lm_step1_mc$r.squared)
  
  # Compute semipartial r between MDIB score and outcome (at Step 2)
  
  ss_total <- anova(lm_step0)$`Sum Sq`
  ss_diff  <- anova_diff$`Sum of Sq`[2]
  
  mdib_scale_sp_rsq_lm_step2 <- ss_diff/ss_total
  mdib_scale_sp_r_lm_step2   <- sqrt(mdib_scale_sp_rsq_lm_step2)
  
  # TODO: Compute semipartial r between BBSIQ score and outcome (at Step 2) given 
  # evidence of mutual suppression (e.g., for incremental predictive validity analysis
  # involving MDIB negative internal bias and NeuroQoL Anxiety). Do this by computing
  # the square root of the R-squared change from Model 1 to Model 2 when Model 1 has
  # only MDIB as a predictor and Model 2 has both MDIB and BBSIQ as predictors (see
  # "2023.11.13 Suppression Notes.docx" for reference).
  
  
  
  
  
  # TODO: Consider computing CIs for semipartial r for MDIB score. However, in
  # simple linear regression, CIs from "r2beta()" (even when "method = 'lm'") differ 
  # from those obtained from "cor.test()", which may use a different method to 
  # compute CIs (even though point estimates are the same). Emailed Byron Jaeger 
  # (developer of "r2glmm" package) about this discrepancy on 11/12/23. Also, in
  # multiple linear regression, point estimates from "r2beta()" (even when "method
  # = 'lm'") differ from those computed manually (not yet sure why).
  
  compute_sp_r(lm_step2_mc)
  
  
  
  
  
  # Collect results in list
  
  res <- list(lm_step0                   = lm_step0,
              lm_step1                   = lm_step1,
              lm_step2                   = lm_step2,
              lm_step1_mc                = lm_step1_mc,
              lm_step2_mc                = lm_step2_mc,
              anova_diff                 = anova_diff,
              summary_lm_step1_mc        = summary_lm_step1_mc,
              summary_lm_step2_mc        = summary_lm_step2_mc,
              confint_lm_step1_mc        = confint_lm_step1_mc,
              confint_lm_step2_mc        = confint_lm_step2_mc,
              bbsiq_scale_r_lm_step1     = bbsiq_scale_r_lm_step1,
              ss_total                   = ss_total,
              ss_diff                    = ss_diff,
              mdib_scale_sp_rsq_lm_step2 = mdib_scale_sp_rsq_lm_step2,
              mdib_scale_sp_r_lm_step2   = mdib_scale_sp_r_lm_step2)
  
  return(res)
}

# Define function to create hierarchical multiple regression table

create_hmr_tbl <- function(res, outcome) {
  df_step1 <- cbind(as.data.frame(res$summary_lm_step1_mc$coefficients),
                    res$confint_lm_step1_mc)
  df_step2 <- cbind(as.data.frame(res$summary_lm_step2_mc$coefficients),
                    res$confint_lm_step2_mc)
  
  df_step1$Outcome <- outcome
  df_step2$Outcome <- outcome
  
  df_step1$Model <- "Step 1"
  df_step2$Model <- "Step 2"
  
  df_step1$Effect <- row.names(df_step1)
  df_step2$Effect <- row.names(df_step2)
  
  target_cols <- c("Outcome", "Model", "Effect")
  df_step1 <- df_step1[, c(target_cols, names(df_step1)[!(names(df_step1) %in% target_cols)])]
  df_step2 <- df_step2[, c(target_cols, names(df_step2)[!(names(df_step2) %in% target_cols)])]
  
  df_step1 <- df_step1[row.names(df_step1) != "(Intercept)", ]
  df_step2 <- df_step2[row.names(df_step2) != "(Intercept)", ]
  
  row.names(df_step1) <- 1:nrow(df_step1)
  row.names(df_step2) <- 1:nrow(df_step2)
  
  df_step1$r <- res$bbsiq_scale_r_lm_step1
  df_step2$r <- NA
  
  df_step1$sp_r <- NA
  
  df_step2$sp_r <- NA
  df_step2$sp_r[grepl("mdib", df_step2$Effect)] <- res$mdib_scale_sp_r_lm_step2
  
  df <- rbind(df_step1, df_step2)
  
  # Rename columns
  
  names(df)[names(df) == "Estimate"]   <- "b"
  names(df)[names(df) == "Std. Error"] <- "SE"
  names(df)[names(df) == "t value"]    <- "t"
  names(df)[names(df) == "Pr(>|t|)"]   <- "p"
  
  # Identify significant parameters
  
  df$sig <- NA
  df$sig <- ifelse(df$`2.5 %` > 0 | df$`97.5 %` < 0, 1, 0)
  
  # Round selected columns and ensure two decimal digits are printed
  
  two_digit_vars   <- c("b", "SE", "t", "2.5 %", "97.5 %", "r", "sp_r")
  three_digit_vars <- "p"
  
  df[, two_digit_vars]   <- format(round(df[, two_digit_vars],   2),
                                   nsmall = 2, trim = TRUE)
  df[, three_digit_vars] <- format(round(df[, three_digit_vars], 3),
                                   nsmall = 3, trim = TRUE)
  
  # Format combined b and SE column
  
  df$`b (SE)` <- paste0(df$b, " (", df$SE, ")")
  
  # Format CIs
  
  df$`95% CI` <- paste0("[", df$`2.5 %`, ", ", df$`97.5 %`, "]")
  
  # Remove redundant columns
  
  rm_cols <- c("b", "SE", "2.5 %", "97.5 %")
  
  df <- df[, !(names(df) %in% rm_cols)]
  
  # Rearrange columns
  
  df <- df[, c("Outcome", "Model", "Effect", "b (SE)", "t", "p", "95% CI", "r", "sp_r", "sig")]
  
  return(df)
}

# ---------------------------------------------------------------------------- #
# Import data ----
# ---------------------------------------------------------------------------- #

load("~/data/further_clean/mdib_hd_dat2.RData")

# ---------------------------------------------------------------------------- #
# Create time point variable, restrict columns, convert to wide format ----
# ---------------------------------------------------------------------------- #

# Create time point variable

mdib_hd_dat2$wave <- NA
mdib_hd_dat2$wave[mdib_hd_dat2$redcap_event_name == "baseline_arm_1"] <- "baseline"
mdib_hd_dat2$wave[mdib_hd_dat2$redcap_event_name == "followup_arm_1"] <- "followup"

# Restrict to columns of interest

anlys_cols <- c("mdib_neg_9_int_m", "mdib_neg_9_ext_m", "bbsiq_neg_int_m", "bbsiq_neg_ext_m", 
                "asi_m", "asi_red_phy_m", "asi_red_cog_m", "asi_red_soc_m", 
                "bfne2_8_m", "neuroqol_anx_m", "sads_m", "sads_red_m", "auditc_m")

target_cols <- c("record_id", "wave", anlys_cols)

mdib_hd_dat2_short <- mdib_hd_dat2[, target_cols]

# Convert to wide format

mdib_hd_dat2_short_wide <- reshape(mdib_hd_dat2_short, direction = "wide",
                                   idvar = "record_id", timevar = "wave",
                                   v.names = anlys_cols)

# Remove columns not assessed at certain time points

rm_cols <- c("sads_red_m.baseline", 
             "asi_m.followup", 
             "asi_red_phy_m.followup", "asi_red_cog_m.followup", "asi_red_soc_m.followup", 
             "bfne2_8_m.followup", "sads_m.followup", "auditc_m.followup")

mdib_hd_dat2_short_wide <- mdib_hd_dat2_short_wide[, !(names(mdib_hd_dat2_short_wide) %in% rm_cols)]

# Export data

save(mdib_hd_dat2_short_wide, file = "~/data/further_clean/mdib_hd_dat2_short_wide.RData")

# ---------------------------------------------------------------------------- #
# Analyze test-retest reliability ----
# ---------------------------------------------------------------------------- #

# Compute 2-week test-retest reliability correlations for MDIB and BBSIQ scales
# using average item scores and pairwise deletion

metric <- "test_retest"

x_var  <- "mdib_neg_9_int_m.baseline"
y_vars <- "mdib_neg_9_int_m.followup"

mdib_neg_9_int_m_retest <- cor_test_one_to_many(mdib_hd_dat2_short_wide, x_var, y_vars, metric)

x_var  <- "mdib_neg_9_ext_m.baseline"
y_vars <- "mdib_neg_9_ext_m.followup"

mdib_neg_9_ext_m_retest <- cor_test_one_to_many(mdib_hd_dat2_short_wide, x_var, y_vars, metric)

x_var  <- "bbsiq_neg_int_m.baseline"
y_vars <- "bbsiq_neg_int_m.followup"

bbsiq_neg_int_m_retest <- cor_test_one_to_many(mdib_hd_dat2_short_wide, x_var, y_vars, metric)

x_var  <- "bbsiq_neg_ext_m.baseline"
y_vars <- "bbsiq_neg_ext_m.followup"

bbsiq_neg_ext_m_retest <- cor_test_one_to_many(mdib_hd_dat2_short_wide, x_var, y_vars, metric)

# Compare dependent, non-overlapping correlations for MDIB versus BBSIQ scales.
# In neither comparison do the correlations significantly differ.

cocor(~mdib_neg_9_int_m.baseline + mdib_neg_9_int_m.followup | 
        bbsiq_neg_int_m.baseline + bbsiq_neg_int_m.followup, 
      mdib_hd_dat2_short_wide, test = "zou2007")

cocor(~mdib_neg_9_ext_m.baseline + mdib_neg_9_ext_m.followup | 
        bbsiq_neg_ext_m.baseline + bbsiq_neg_ext_m.followup, 
      mdib_hd_dat2_short_wide, test = "zou2007")

# ---------------------------------------------------------------------------- #
# Analyze convergent validity ----
# ---------------------------------------------------------------------------- #

# Analyze correlations at baseline using average item scores and pairwise deletion

metric <- "convergent"

x_var  <- "mdib_neg_9_int_m.baseline"
y_vars <- c("bbsiq_neg_int_m.baseline",
            "asi_m.baseline", "asi_red_phy_m.baseline", "asi_red_cog_m.baseline")

mdib_neg_9_int_m_conv <- cor_test_one_to_many(mdib_hd_dat2_short_wide, x_var, y_vars, metric)

x_var  <- "mdib_neg_9_ext_m.baseline"
y_vars <- c("bbsiq_neg_ext_m.baseline", 
            "asi_m.baseline", "asi_red_soc_m.baseline",
            "bfne2_8_m.baseline")

mdib_neg_9_ext_m_conv <- cor_test_one_to_many(mdib_hd_dat2_short_wide, x_var, y_vars, metric)

# ---------------------------------------------------------------------------- #
# Analyze discriminant validity ----
# ---------------------------------------------------------------------------- #

# Analyze correlations at baseline using average item scores and pairwise deletion

metric <- "discriminant"
y_vars <- "auditc_m.baseline"

x_var  <- "mdib_neg_9_int_m.baseline"

mdib_neg_9_int_m_disc <- cor_test_one_to_many(mdib_hd_dat2_short_wide, x_var, y_vars, metric)

x_var  <- "mdib_neg_9_ext_m.baseline"

mdib_neg_9_ext_m_disc <- cor_test_one_to_many(mdib_hd_dat2_short_wide, x_var, y_vars, metric)

# ---------------------------------------------------------------------------- #
# Analyze concurrent validity ----
# ---------------------------------------------------------------------------- #

# Analyze correlations at baseline using average item scores and pairwise deletion

metric <- "concurrent"

x_var  <- "mdib_neg_9_int_m.baseline"
y_vars <- "neuroqol_anx_m.baseline"

mdib_neg_9_int_m_conc <- cor_test_one_to_many(mdib_hd_dat2_short_wide, x_var, y_vars, metric)

x_var  <- "mdib_neg_9_ext_m.baseline"
y_vars <- c("neuroqol_anx_m.baseline", "sads_m.baseline")

mdib_neg_9_ext_m_conc <- cor_test_one_to_many(mdib_hd_dat2_short_wide, x_var, y_vars, metric)

# ---------------------------------------------------------------------------- #
# Analyze predictive validity ----
# ---------------------------------------------------------------------------- #

# Analyze correlations from baseline to follow-up using average item scores
# and pairwise deletion

metric <- "predictive"

x_var  <- "mdib_neg_9_int_m.baseline"
y_vars <- "neuroqol_anx_m.followup"

mdib_neg_9_int_m_pred <- cor_test_one_to_many(mdib_hd_dat2_short_wide, x_var, y_vars, metric)

x_var  <- "mdib_neg_9_ext_m.baseline"
y_vars <- c("neuroqol_anx_m.followup", "sads_red_m.followup")

mdib_neg_9_ext_m_pred <- cor_test_one_to_many(mdib_hd_dat2_short_wide, x_var, y_vars, metric)

# ---------------------------------------------------------------------------- #
# Create table for convergent, discriminant, concurrent, and predictive validity ----
# ---------------------------------------------------------------------------- #

tbl_external_validity <- rbind(mdib_neg_9_int_m_conv, mdib_neg_9_ext_m_conv, 
                               mdib_neg_9_int_m_disc, mdib_neg_9_ext_m_disc, 
                               mdib_neg_9_int_m_conc, mdib_neg_9_ext_m_conc, 
                               mdib_neg_9_int_m_pred, mdib_neg_9_ext_m_pred)

# ---------------------------------------------------------------------------- #
# Analyze incremental concurrent validity ----
# ---------------------------------------------------------------------------- #

# Given no missingness in predictors, each analysis is restricted to participants with
# outcome (i.e., same sample size as concurrent validity analyses above)

predictors <- c("bbsiq_neg_int_m.baseline", "bbsiq_neg_ext_m.baseline",
                "mdib_neg_9_int_m.baseline", "mdib_neg_9_ext_m.baseline")

sum(is.na(mdib_hd_dat2_short_wide[predictors])) == 0

# Run function to analyze incremental concurrent validity

res_conc_neuroqol_mdib_neg_9_int <- run_hmr(mdib_hd_dat2_short_wide,
                                            "neuroqol_anx_m.baseline",
                                            "bbsiq_neg_int_m.baseline", "mdib_neg_9_int_m.baseline")

res_conc_neuroqol_mdib_neg_9_ext <- run_hmr(mdib_hd_dat2_short_wide,
                                            "neuroqol_anx_m.baseline",
                                            "bbsiq_neg_ext_m.baseline", "mdib_neg_9_ext_m.baseline")
res_conc_sads_mdib_neg_9_ext     <- run_hmr(mdib_hd_dat2_short_wide,
                                            "sads_m.baseline",
                                            "bbsiq_neg_ext_m.baseline", "mdib_neg_9_ext_m.baseline")

# Create tables

tbl_conc_neuroqol_mdib_neg_9_int <- create_hmr_tbl(res_conc_neuroqol_mdib_neg_9_int, "neuroqol_anx_m.baseline")

tbl_conc_neuroqol_mdib_neg_9_ext <- create_hmr_tbl(res_conc_neuroqol_mdib_neg_9_ext, "neuroqol_anx_m.baseline")
tbl_conc_sads_mdib_neg_9_ext     <- create_hmr_tbl(res_conc_sads_mdib_neg_9_ext,     "sads_m.baseline")

tbl_conc <- rbind(tbl_conc_neuroqol_mdib_neg_9_int, 
                  tbl_conc_neuroqol_mdib_neg_9_ext, tbl_conc_sads_mdib_neg_9_ext)

tbl_conc$metric <- "incremental concurrent"
tbl_conc <- tbl_conc[, c("metric", names(tbl_conc)[!(names(tbl_conc) %in% "metric")])]

# ---------------------------------------------------------------------------- #
# Analyze incremental predictive validity ----
# ---------------------------------------------------------------------------- #

# Given no missingness in predictors (see above), each analysis is restricted to 
# participants with outcome (same sample size as predictive validity analyses above)

# Run function to analyze incremental predictive validity

res_pred_neuroqol_mdib_neg_9_int <- run_hmr(mdib_hd_dat2_short_wide,
                                            "neuroqol_anx_m.followup",
                                            "bbsiq_neg_int_m.baseline", "mdib_neg_9_int_m.baseline")

res_pred_neuroqol_mdib_neg_9_ext <- run_hmr(mdib_hd_dat2_short_wide,
                                            "neuroqol_anx_m.followup",
                                            "bbsiq_neg_ext_m.baseline", "mdib_neg_9_ext_m.baseline")
res_pred_sads_mdib_neg_9_ext     <- run_hmr(mdib_hd_dat2_short_wide,
                                            "sads_red_m.followup",
                                            "bbsiq_neg_ext_m.baseline", "mdib_neg_9_ext_m.baseline")

# Create tables

tbl_pred_neuroqol_mdib_neg_9_int <- create_hmr_tbl(res_pred_neuroqol_mdib_neg_9_int, "neuroqol_anx_m.followup")

tbl_pred_neuroqol_mdib_neg_9_ext <- create_hmr_tbl(res_pred_neuroqol_mdib_neg_9_ext, "neuroqol_anx_m.followup")
tbl_pred_sads_mdib_neg_9_ext     <- create_hmr_tbl(res_pred_sads_mdib_neg_9_ext,     "sads_red_m.followup")

tbl_pred <- rbind(tbl_pred_neuroqol_mdib_neg_9_int, 
                  tbl_pred_neuroqol_mdib_neg_9_ext, tbl_pred_sads_mdib_neg_9_ext)

tbl_pred$metric <- "incremental predictive"
tbl_pred <- tbl_pred[, c("metric", names(tbl_pred)[!(names(tbl_pred) %in% "metric")])]

# ---------------------------------------------------------------------------- #
# Create table for incremental validity ----
# ---------------------------------------------------------------------------- #

tbl_incremental_validity <- rbind(tbl_conc, tbl_pred)

# ---------------------------------------------------------------------------- #
# Set "flextable" defaults and section and text properties ----
# ---------------------------------------------------------------------------- #

# "flextable" defaults are set in "set_flextable_defaults.R" above

# Section and text properties are sourced from "set_officer_properties.R" above

# ---------------------------------------------------------------------------- #
# Format table for convergent, discriminant, concurrent, and predictive validity ----
# ---------------------------------------------------------------------------- #

# "flextable" defaults are set in "set_flextable_defaults.R" above

# Define function to format external validity table (except incremental validity)

format_tbl_external_validity <- function(tbl_external_validity, gen_note, title) {
  tbl <- tbl_external_validity
  
  left_align_body_cols <- c("x_var", "y_var")
  
  target_cols <- c("metric", "x_var", "y_var", "r", "t", "df", "p", "ci", "n")

  # Restructure data as grouped data for table spanners based on "metric"
  
  tbl <- as_grouped_data(tbl[, target_cols], groups = "metric",
                         columns = target_cols[!(target_cols %in% "metric")])
  
  # Define visible columns
  
  visible_cols <- target_cols[!(target_cols %in% "metric")]

  # Create flextable
  
  tbl_ft <- as_flextable(tbl, col_keys = visible_cols) |>
    set_table_properties(align = "left") |>
    
    set_caption(as_paragraph(as_i(title)), word_stylename = "heading 1",
                fp_p = fp_par(padding.left = 0, padding.right = 0),
                align_with_table = FALSE) |>
    
    align(align = "center", part = "header") |>
    align(align = "center", part = "body") |>
    align(j = left_align_body_cols, align = "left", part = "body") |>
    align(align = "left", part = "footer") |>
    
    align(j = 1, i = ~ !is.na(metric), align = "center", part = "body") |>
    hline(j = 1, i = ~ !is.na(metric), part = "body") |>
    bold(j = 1,  i = ~ !is.na(metric), part = "body") |>
    
    set_header_labels(x_var = "MDIB Measure",
                      y_var = "Validity Measure",
                      ci    = "95% CI") |>
    compose(j = "r",  part = "header", value = as_paragraph(as_i("r"))) |>
    compose(j = "t",  part = "header", value = as_paragraph(as_i("t"))) |>
    compose(j = "df", part = "header", value = as_paragraph(as_i("df"))) |>
    compose(j = "p",  part = "header", value = as_paragraph(as_i("p"))) |>
    compose(j = "n",  part = "header", value = as_paragraph(as_i("n"))) |>
    
    labelizor(part = "body",
              labels = c("mdib_neg_9_int_m.baseline" = "Negative Internal Bias",
                         "mdib_neg_9_ext_m.baseline" = "Negative External Bias",
                         "bbsiq_neg_int_m.baseline"  = "Negative Internal Bias (BBSIQ)",
                         "bbsiq_neg_ext_m.baseline"  = "Negative External Bias (BBSIQ)",
                         "asi_m.baseline"            = "Overall Anxiety Sensitivity (ASI)",
                         "asi_red_phy_m.baseline"    = "Anxiety Sensitivity-Physical (2 ASI items)",
                         "asi_red_cog_m.baseline"    = "Anxiety Sensitivity-Cognitive (2 ASI items)",
                         "asi_red_soc_m.baseline"    = "Anxiety Sensitivity-Social (2 ASI items)",
                         "bfne2_8_m.baseline"        = "Fear of Negative Evaluation (8 BFNE-II items)",
                         "auditc_m.baseline"         = "Alcohol Use (AUDIT-C)",
                         "neuroqol_anx_m.baseline"   = "Anxiety Symptoms (Neuro-QoL Anxiety)",
                         "neuroqol_anx_m.followup"   = "Anxiety Symptoms (Neuro-QoL Anxiety)",
                         "sads_m.baseline"           = "Social Avoidance and Distress (SADS)",
                         "sads_red_m.followup"       = "Social Avoidance and Distress (8 SADS items)")) |>

    add_footer_lines(gen_note) |>
    autofit()
  
  return(tbl_ft)
}

# Define general note

gen_note <- as_paragraph_md("*Note.* MDIB measures (3 negative ratings for internal threats; 6 negative ratings for external threats) and convergent, discriminant, and concurrent validity measures were assessed at baseline; predictive validity measures were assessed at 2-week follow-up. MDIB = Movement Disorders Interpretation Bias Scale; BBSIQ = Brief Body Sensations Interpretation Questionnaire; ASI = Anxiety Sensitivity Index; BFNE-II = Brief Fear of Negative Evaluation-II Scale; AUDIT-C = Alcohol Use Disorders Identification Test-Consumption Items; Neuro-QoL Anxiety = Neuro-QoL Short Form for Anxiety; SADS = Social Avoidance and Distress Scale.")

# Run function

tbl_external_validity_ft <- format_tbl_external_validity(tbl_external_validity, gen_note,
  "Associations of MDIB Negative Bias Measures With Convergent, Discriminant, Concurrent, and Predictive Validity Measures")

# ---------------------------------------------------------------------------- #
# Format table for incremental validity ----
# ---------------------------------------------------------------------------- #

# Define function to format incremental validity table

format_tbl_incremental_validity <- function(tbl_incremental_validity, gen_note, title) {
  tbl <- tbl_incremental_validity
  
  left_align_body_cols <- c("Outcome", "Model", "Effect")
  
  target_cols <- c("metric", "Outcome", "Model", "Effect", "b (SE)", "t", "p", 
                   "95% CI", "r", "sp_r", "sig")
  
  # Restructure data as grouped data for table spanners based on "metric"
  
  tbl <- as_grouped_data(tbl[, target_cols], groups = "metric",
                         columns = target_cols[!(target_cols %in% "metric")])
  
  # Define visible columns
  
  visible_cols <- target_cols[!(target_cols %in% c("metric", "sig"))]
  
  # Create flextable
  
  tbl_ft <- as_flextable(tbl, col_keys = visible_cols) |>
    set_table_properties(align = "left") |>
    
    set_caption(as_paragraph(as_i(title)), word_stylename = "heading 1",
                fp_p = fp_par(padding.left = 0, padding.right = 0),
                align_with_table = FALSE) |>
    
    align(align = "center", part = "header") |>
    align(align = "center", part = "body") |>
    align(j = left_align_body_cols, align = "left", part = "body") |>
    align(align = "left", part = "footer") |>
    
    align(j = 1, i = ~ !is.na(metric), align = "center", part = "body") |>
    hline(j = 1, i = ~ !is.na(metric), part = "body") |>
    bold(j = 1,  i = ~ !is.na(metric), part = "body") |>
    
    set_header_labels(x_var = "MDIB Measure",
                      y_var = "Validity Measure",
                      ci    = "95% CI") |>
    compose(j = "b (SE)", part = "header", 
            value = as_paragraph(as_i("b"), " (", as_i("SE"), ")")) |>
    compose(j = "t",    part = "header", value = as_paragraph(as_i("t"))) |>
    compose(j = "p",    part = "header", value = as_paragraph(as_i("p"))) |>
    compose(j = "r",    part = "header", value = as_paragraph(as_i("r"))) |>
    compose(j = "sp_r", part = "header", value = as_paragraph(as_i("r"), as_sub("sp"))) |>
    
    labelizor(part = "body",
              labels = c("mdib_neg_9_int_m.baselinec" = "Neg. Int. Bias (MDIB)",
                         "mdib_neg_9_ext_m.baselinec" = "Neg. Ext. Bias (MDIB)",
                         "bbsiq_neg_int_m.baselinec"  = "Neg. Int. Bias (BBSIQ)",
                         "bbsiq_neg_ext_m.baselinec"  = "Neg. Ext. Bias (BBSIQ)",
                         "neuroqol_anx_m.baseline"   = "Anxiety Symptoms (Neuro-QoL Anxiety)",
                         "neuroqol_anx_m.followup"   = "Anxiety Symptoms (Neuro-QoL Anxiety)",
                         "sads_m.baseline"           = "Social Avoid./Distress (SADS)",
                         "sads_red_m.followup"       = "Social Avoid./Distress (8 SADS items)",
                         "NA"                        = "")) |>
    
    add_footer_lines(gen_note) |>
    autofit()
  
  return(tbl_ft)
}

# Define general note

gen_note <- as_paragraph_md("*Note.* Key rows testing the incremental validity of the MDIB are in boldface. MDIB measures (3 negative ratings for internal threats; 6 negative ratings for external threats), BBSIQ measures, and concurrent validity outcomes were assessed at baseline; predictive validity outcomes were assessed at 2-week follow-up. MDIB = Movement Disorders Interpretation Bias Scale; BBSIQ = Brief Body Sensations Interpretation Questionnaire; Neuro-QoL Anxiety = Neuro-QoL Short Form for Anxiety; SADS = Social Avoidance and Distress Scale.")

# Run function

tbl_incremental_validity_ft <- 
  format_tbl_incremental_validity(tbl_incremental_validity, gen_note,
  "Hierarchical Multiple Regression Models Testing Incremental Concurrent and Predictive Validity of MDIB Negative Bias Measures")

# ---------------------------------------------------------------------------- #
# Write external validity tables to MS Word ----
# ---------------------------------------------------------------------------- #

# Write external validity tables (note: "flextable" seems to have a bug in which 
# blank page is at end of doc)

external_validity_tbls <- list(tbl_external_validity_ft, tbl_incremental_validity_ft)

external_validity_tbls_orientations <- c("l", "l")
external_validity_tbls_numbers      <- c("3", "4")

doc <- read_docx()
doc <- body_set_default_section(doc, psect_prop)

for (i in 1:length(external_validity_tbls)) {
  doc <- body_add_fpar(doc, fpar(ftext(paste0("Table ", external_validity_tbls_numbers[[i]]),
                                       prop = text_prop_bold)))
  doc <- body_add_par(doc, "")
  
  doc <- body_add_flextable(doc, external_validity_tbls[[i]], align = "left")
  
  if (external_validity_tbls_orientations[[i]] == "p") {
    doc <- body_end_block_section(doc, block_section(psect_prop))
  } else if (external_validity_tbls_orientations[[i]] == "l") {
    doc <- body_end_block_section(doc, block_section(lsect_prop))
  }
}

external_validity_tbl_path <- "./results/external_validity/tables/"

dir.create(external_validity_tbl_path, recursive = TRUE)

print(doc, target = paste0(external_validity_tbl_path, "external_validity_tbls.docx"))
