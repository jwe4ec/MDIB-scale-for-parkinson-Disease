# ---------------------------------------------------------------------------- #
# Compute Internal Consistency
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

source("./code/1a_define_functions.R")

# Check correct R version, load groundhog package, and specify groundhog_day

groundhog_day <- version_control()

# Load packages

pkgs <- c("lavaan", "semTools", "MBESS")
groundhog.library(pkgs, groundhog_day)

# Set seed

set.seed(1234)

# ---------------------------------------------------------------------------- #
# Define functions used in script ----
# ---------------------------------------------------------------------------- #

# Define function to export internal consistency results (except CIs)

export_omega_res <- function(summary, residuals, reliability, omega, 
                             path, filename) {
  sink(paste0(path, paste0(filename, ".txt")))
  cat("Summary:", "\n", "\n")
  print(summary)
  cat("Residual Correlations:", "\n", "\n")
  print(residuals)
  cat("Output from reliability():", "\n", "\n")
  print(reliability)
  cat("\n")
  cat("Output from compRelSEM():", "\n", "\n")
  print(omega)
  sink()
}

# Define function to export internal consistency CIs

export_omega_ci_res <- function(omega_ci, path, filename) {
  sink(paste0(path, paste0(filename, ".txt")))
  cat("Output from ci.reliability():", "\n", "\n")
  print(omega_ci)
  sink()
}

# ---------------------------------------------------------------------------- #
# Import data ----
# ---------------------------------------------------------------------------- #

load("~/data/further_clean/mdib_bl_neg_9.Rdata")
load("~/data/further_clean/mdib_bl_neg_9_ord.Rdata")

# ---------------------------------------------------------------------------- #
# Define internal and external scales for 9 reduced negative bias items ----
# ---------------------------------------------------------------------------- #

mdib_neg_9_int_items <- c("mdib_neg_int_cleaning_4c", 
                          "mdib_neg_int_medication_8b", "mdib_neg_int_cough_12b")

mdib_neg_9_ext_items <- c("mdib_neg_ext_reminder_3c", "mdib_neg_ext_neighbor_5a", 
                          "mdib_neg_ext_exercise_7a", "mdib_neg_ext_walk_9c",
                          "mdib_neg_ext_job_10a", "mdib_neg_ext_stumble_11c")

mdib_neg_9_items <- list(mdib_neg_9_int_items = mdib_neg_9_int_items,
                         mdib_neg_9_ext_items = mdib_neg_9_ext_items)

# ---------------------------------------------------------------------------- #
# Compute internal consistency for negative internal items ----
# ---------------------------------------------------------------------------- #

# Prepare lavaan syntax for latent factor

items_fml <- paste(mdib_neg_9_items$mdib_neg_9_int_items, collapse = " + ")

# Compute omega for negative internal items at baseline based on single-factor 
# model following Flora (2020). Use "MLR" estimator per our preregistration and
# Flora (2020), but also use "MLM" and "WLSMV" as these were used for the EFA

mod_mdib_neg_9_int <-
  'f1 =~  mdib_neg_int_cleaning_4c + 
         mdib_neg_int_medication_8b + mdib_neg_int_cough_12b'

fit_mdib_neg_9_int_mlr   <- cfa(mod_mdib_neg_9_int, 
                                data = mdib_bl_neg_9,
                                std.lv = TRUE,
                                estimator = "MLR")

fit_mdib_neg_9_int_mlm   <- cfa(mod_mdib_neg_9_int, 
                                data = mdib_bl_neg_9,
                                std.lv = TRUE,
                                estimator = "MLM")

fit_mdib_neg_9_int_wlsmv <- cfa(mod_mdib_neg_9_int, 
                                data = mdib_bl_neg_9_ord,
                                std.lv = TRUE,
                                ordered = TRUE,
                                estimator = "WLSMV")

# Interpreting CFA fit based on Flora (2020), robust CFI is >= .90, robust TLI is 
# not >= 0.90, and robust RMSEA is not <= .08. 
# DANDAN: all NA, CFI = TLI=1, RMSEA=0 

summ_fit_mdib_neg_9_int_mlr   <- summary(fit_mdib_neg_9_int_mlr, fit.measures = TRUE)
summ_fit_mdib_neg_9_int_mlm   <- summary(fit_mdib_neg_9_int_mlm, fit.measures = TRUE)
summ_fit_mdib_neg_9_int_wlsmv <- summary(fit_mdib_neg_9_int_wlsmv, fit.measures = TRUE)



resid_fit_mdib_neg_9_int_mlr   <- residuals(fit_mdib_neg_9_int_mlr, type = 'cor')
resid_fit_mdib_neg_9_int_mlm   <- residuals(fit_mdib_neg_9_int_mlm, type = 'cor')
resid_fit_mdib_neg_9_int_wlsmv <- residuals(fit_mdib_neg_9_int_wlsmv, type = 'cor')

# Note: "reliability" function is now deprecated. "compRelSEM" yields identical
# results as "omega2" (based on model-implied variance) from "reliability".
#   For MLR,   .794 (note: this matches point estimate from "ci.reliability" below,
#                    whereas "omega3", based on observed sample variance, is .802)
#   For MLM,   .794 (note: "omega3" is .802)
#   For WLSMV, .804 (note: "omega3" is .816, which in this case is what matches 
#                    point estimate from "ci.reliability" below)

#reliability_fit_mdib_neg_9_int_mlr   <- reliability(fit_mdib_neg_9_int_mlr)
#reliability_fit_mdib_neg_9_int_mlm   <- reliability(fit_mdib_neg_9_int_mlm)
#reliability_fit_mdib_neg_9_int_wlsmv <- reliability(fit_mdib_neg_9_int_wlsmv)
#In reliability(fit_mdib_neg_9_int_wlsmv) : 
#The reliability() function was deprecated in 2022 and will cease to be included in future versions of semTools. See help('semTools-deprecated) for details.

#It is replaced by the compRelSEM() function, which can estimate alpha and model-based reliability in an even wider variety of models and data types, with greater control in specifying the desired type of reliability coefficient (i.e., more explicitly choosing assumptions). 

#The average variance extracted should never have been included because it is not a reliability coefficient. It is now available from the AVE() function.
omega_fit_mdib_neg_9_int_mlr <- 
  compRelSEM(fit_mdib_neg_9_int_mlr,   obs.var = FALSE, tau.eq = FALSE, ord.scale = FALSE)
omega_fit_mdib_neg_9_int_mlm <- 
  compRelSEM(fit_mdib_neg_9_int_mlm,   obs.var = FALSE, tau.eq = FALSE, ord.scale = FALSE)
omega_fit_mdib_neg_9_int_wlsmv <- 
  compRelSEM(fit_mdib_neg_9_int_wlsmv, obs.var = FALSE, tau.eq = FALSE, ord.scale = TRUE)

# Computing 95% percentile bootstrap CIs yields warnings of this type:
#   In lav_object_post_check(object) :
#     lavaan WARNING: some estimated ov variances are negative

# Results:
#   omega             = .761, 
#   categorical omega = .76, 

mdib_bl_neg_9_int <- mdib_bl_neg_9[, mdib_neg_9_int_items]
mdib_bl_neg_9_int_ord <- mdib_bl_neg_9_ord[, mdib_neg_9_int_items]

# Note: Same results to 2 decimal places with more (5000) iterations (see Dunn et 
# al., 2014, p. 408) for omega, so use 1000 for omega_cat and other omegas below

omega_ci_mdib_bl_neg_9_int_1000iter <- ci.reliability(data = mdib_bl_neg_9_int,
                                                      type = "omega",
                                                      conf.level = 0.95,
                                                      interval.type = "perc",
                                                      B = 1000)

omega_ci_mdib_bl_neg_9_int_ord_1000iter <- ci.reliability(data = mdib_bl_neg_9_int_ord,
                                                          type = "categorical",
                                                          conf.level = 0.95,
                                                          interval.type = "perc",
                                                          B = 1000)

# Export results (except CIs)

internal_consistency_path   <- "./results/internal_consistency/"

dir.create(internal_consistency_path, recursive = TRUE)

export_omega_res(summ_fit_mdib_neg_9_int_mlr, resid_fit_mdib_neg_9_int_mlr,
                 reliability_fit_mdib_neg_9_int_mlr, omega_fit_mdib_neg_9_int_mlr,
                 internal_consistency_path, "int_omega_cfa_mlr")
export_omega_res(summ_fit_mdib_neg_9_int_mlm, resid_fit_mdib_neg_9_int_mlm,
                 reliability_fit_mdib_neg_9_int_mlm, omega_fit_mdib_neg_9_int_mlm,
                 internal_consistency_path, "int_omega_cfa_mlm")
export_omega_res(summ_fit_mdib_neg_9_int_wlsmv, resid_fit_mdib_neg_9_int_wlsmv,
                 reliability_fit_mdib_neg_9_int_wlsmv, omega_fit_mdib_neg_9_int_wlsmv,
                 internal_consistency_path, "int_omega_cfa_wlsmv")

# Export CIs

export_omega_ci_res(omega_ci_mdib_bl_neg_9_int_1000iter,
                    internal_consistency_path, "int_omega_ci")
export_omega_ci_res(omega_ci_mdib_bl_neg_9_int_ord_1000iter,
                    internal_consistency_path, "int_omega_ci_cat")

# ---------------------------------------------------------------------------- #
# Compute internal consistency for negative external items ----
# ---------------------------------------------------------------------------- #

# Prepare lavaan syntax for latent factor

items_fml <- paste(mdib_neg_9_items$mdib_neg_9_ext_items, collapse = " + ")

# Compute omega for negative external items at baseline based on single-factor 
# model following Flora (2020). Use "MLR" estimator per our preregistration and
# Flora (2020), but also use "MLM" and "WLSMV" as these were used for the EFA

mod_mdib_neg_9_ext <-
  'f1 =~ mdib_neg_ext_reminder_3c + mdib_neg_ext_neighbor_5a + 
         mdib_neg_ext_exercise_7a + mdib_neg_ext_walk_9c + mdib_neg_ext_job_10a + mdib_neg_ext_stumble_11c'

fit_mdib_neg_9_ext_mlr   <- cfa(mod_mdib_neg_9_ext, 
                                data = mdib_bl_neg_9,
                                std.lv = TRUE,
                                estimator = "MLR")

fit_mdib_neg_9_ext_mlm   <- cfa(mod_mdib_neg_9_ext, 
                                data = mdib_bl_neg_9,
                                std.lv = TRUE,
                                estimator = "MLM")

fit_mdib_neg_9_ext_wlsmv <- cfa(mod_mdib_neg_9_ext, 
                                data = mdib_bl_neg_9_ord,
                                std.lv = TRUE,
                                ordered = TRUE,
                                estimator = "WLSMV")

# Interpreting CFA fit based on Flora (2020), robust CFI and robust TLI are >= .90
# and robust RMSEA is <= .08

summ_fit_mdib_neg_9_ext_mlr   <- summary(fit_mdib_neg_9_ext_mlr, fit.measures = TRUE)
summ_fit_mdib_neg_9_ext_mlm   <- summary(fit_mdib_neg_9_ext_mlm, fit.measures = TRUE)
summ_fit_mdib_neg_9_ext_wlsmv <- summary(fit_mdib_neg_9_ext_wlsmv, fit.measures = TRUE)


resid_fit_mdib_neg_9_ext_mlr   <- residuals(fit_mdib_neg_9_ext_mlr, type = 'cor')
resid_fit_mdib_neg_9_ext_mlm   <- residuals(fit_mdib_neg_9_ext_mlm, type = 'cor')
resid_fit_mdib_neg_9_ext_wlsmv <- residuals(fit_mdib_neg_9_ext_wlsmv, type = 'cor')

#reliability_fit_mdib_neg_9_ext_mlr   <- reliability(fit_mdib_neg_9_ext_mlr)
reliability_fit_mdib_neg_9_ext_mlm   <- reliability(fit_mdib_neg_9_ext_mlm)
reliability_fit_mdib_neg_9_ext_wlsmv <- reliability(fit_mdib_neg_9_ext_wlsmv)

omega_fit_mdib_neg_9_ext_mlr <-
  compRelSEM(fit_mdib_neg_9_ext_mlr,   obs.var = FALSE, tau.eq = FALSE, ord.scale = FALSE)
omega_fit_mdib_neg_9_ext_mlm <-
  compRelSEM(fit_mdib_neg_9_ext_mlm,   obs.var = FALSE, tau.eq = FALSE, ord.scale = FALSE)
omega_fit_mdib_neg_9_ext_wlsmv <-
  compRelSEM(fit_mdib_neg_9_ext_wlsmv, obs.var = FALSE, tau.eq = FALSE, ord.scale = TRUE)


mdib_bl_neg_9_ext <- mdib_bl_neg_9[, mdib_neg_9_ext_items]
mdib_bl_neg_9_ext_ord <- mdib_bl_neg_9_ord[, mdib_neg_9_ext_items]

omega_ci_mdib_bl_neg_9_ext_1000iter <- ci.reliability(data = mdib_bl_neg_9_ext,
                                                      type = "omega",
                                                      conf.level = 0.95,
                                                      interval.type = "perc",
                                                      B = 1000)

omega_ci_mdib_bl_neg_9_ext_ord_1000iter <- ci.reliability(data = mdib_bl_neg_9_ext_ord,
                                                          type = "categorical",
                                                          conf.level = 0.95,
                                                          interval.type = "perc",
                                                          B = 1000)

# Export results (except CIs)

export_omega_res(summ_fit_mdib_neg_9_ext_mlr, resid_fit_mdib_neg_9_ext_mlr,
                 reliability_fit_mdib_neg_9_ext_mlr, omega_fit_mdib_neg_9_ext_mlr,
                 internal_consistency_path, "ext_omega_cfa_mlr")
export_omega_res(summ_fit_mdib_neg_9_ext_mlm, resid_fit_mdib_neg_9_ext_mlm,
                 reliability_fit_mdib_neg_9_ext_mlm, omega_fit_mdib_neg_9_ext_mlm,
                 internal_consistency_path, "ext_omega_cfa_mlm")
export_omega_res(summ_fit_mdib_neg_9_ext_wlsmv, resid_fit_mdib_neg_9_ext_wlsmv,
                 reliability_fit_mdib_neg_9_ext_wlsmv, omega_fit_mdib_neg_9_ext_wlsmv,
                 internal_consistency_path, "ext_omega_cfa_wlsmv")

# Export CIs

export_omega_ci_res(omega_ci_mdib_bl_neg_9_ext_1000iter,
                    internal_consistency_path, "ext_omega_ci")
export_omega_ci_res(omega_ci_mdib_bl_neg_9_ext_ord_1000iter,
                    internal_consistency_path, "ext_omega_ci_cat")
