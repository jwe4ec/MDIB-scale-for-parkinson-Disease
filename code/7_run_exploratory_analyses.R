# ---------------------------------------------------------------------------- #
# Run Exploratory Analyses
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

# No packages loaded

# Set seed

set.seed(1234)

# ---------------------------------------------------------------------------- #
# Import data ----
# ---------------------------------------------------------------------------- #

load("./data/further_clean/mdib_hd_dat2_short_wide.RData")

# ---------------------------------------------------------------------------- #
# Inspect distribution of Neuro-QoL Anxiety ----
# ---------------------------------------------------------------------------- #

hist(mdib_hd_dat2_short_wide$neuroqol_anx_m.baseline)
hist(mdib_hd_dat2_short_wide$neuroqol_anx_m.followup)

# ---------------------------------------------------------------------------- #
# Compute NeuroQoL-Anxiety total score ----
# ---------------------------------------------------------------------------- #

mdib_hd_dat2_short_wide$neuroqol_anx_tot.baseline <- mdib_hd_dat2_short_wide$neuroqol_anx_m.baseline*8
mdib_hd_dat2_short_wide$neuroqol_anx_tot.followup <- mdib_hd_dat2_short_wide$neuroqol_anx_m.followup*8


x <- as.numeric(as.matrix(mdib_hd_dat2_short_wide$neuroqol_anx_tot.baseline))
mean_all <- mean(x, na.rm = TRUE)
sd_all   <- sd(x, na.rm = TRUE)

mean_all
sd_all


# ---------------------------------------------------------------------------- #
# Compute selected descriptives ----
# ---------------------------------------------------------------------------- #

# Define function for computing sample size, summary, and standard deviation

compute_desc <- function(df, var) {
  n <- sum(!is.na(df[, var]))
  summary <- round(summary(df[, var]), 2)
  sd <- format(round(sd(df[, var]), 2),
               nsmall = 2, trim = TRUE)
  
  cat("Variable:", var, "\n\n")
  
  cat("Summary: ", "\n")
  print(summary)
  cat("\n")
  
  cat("SD =", sd, "\n")
  cat("n =", n, "\n\n")
}

# Run function and export results

dir.create("./results/descriptives")

sink(file = "./results/descriptives/descriptives.txt")
compute_desc(mdib_hd_dat2_short_wide, "neuroqol_anx_m.baseline")
compute_desc(mdib_hd_dat2_short_wide, "neuroqol_anx_tot.baseline")
compute_desc(mdib_hd_dat2_short_wide, "mdib_neg_9_int_m.baseline")
compute_desc(mdib_hd_dat2_short_wide, "mdib_neg_9_ext_m.baseline")
sink()

