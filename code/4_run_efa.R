# ---------------------------------------------------------------------------- #
# Run Exploratory Factor Analysis for the PD paper
# ---------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------- #
# Notes ----
# ---------------------------------------------------------------------------- #

# Before running this script:
# 1. Restart R.
# 2. Set the working directory to the project root folder.
# 3. Make sure the PD data preparation scripts export PD-labeled objects.
#
# This script follows the PD preregistered EFA plan:
# - Inspect item distributions first.
# - If items are severely skewed, treat items as categorical.
# - Use parallel analysis directly, without scree plots.
# - For categorical items:
#   - Use polychoric correlations in parallel analyses.
#   - Use WLSMV estimation for EFAs.
# - Use oblimin rotation as the primary solution.
# - Use geomin and promax rotations as sensitivity analyses.
#
# This script intentionally removes:
# - scree plot analyses,
# - parallel analyses based on Pearson correlations,
# - MLM EFA models,
# - HD-specific comments and decision logic.
#
# Output folders and file stems are intentionally abbreviated to avoid creating
# very long file paths in the results directory.
#
# Result-informed notes from the current PD run:
# - Visual inspection showed substantial skewness, especially for the 12 negative
#   items, which showed pronounced floor effects and sparse upper-category
#   responses. This supported treating the items as ordered categorical indicators.
# - For all 36 items, polychoric parallel analysis suggested an upper bound of
#   three components, so 2-, 3-, and 4-factor WLSMV EFA solutions were examined.
# - The 36-item EFAs did not yield a clear, stable, theoretically interpretable
#   structure beyond a broad benign-versus-negative distinction. Therefore, the
#   analyses proceeded to the 12 theorized negative items, consistent with the
#   analysis plan and the theorized negative-bias item set.
# - For the 12 negative items, polychoric parallel analysis suggested a one-
#   component upper bound. Following the preregistered +/- 1 approach, 1- and
#   2-factor WLSMV EFA solutions were examined.
# - The 12-item 1-factor solution supported a broad negative bias factor but had
#   poor model fit and a very weak item, mdib_neg_int_remember_1b.
# - The 12-item 2-factor solution was more consistent with the theorized
#   internal/external distinction, but three item-level concerns remained:
#     1. mdib_neg_int_remember_1b did not load saliently on either factor.
#     2. mdib_neg_ext_server_2a loaded with the internal rather than external items.
#     3. mdib_neg_int_email_6b showed a mild cross-loading in some rotations.
# - Across the item-removal sequences, mdib_neg_int_remember_1b and
#   mdib_neg_ext_server_2a were the clearest problematic items. The evidence for
#   removing mdib_neg_int_email_6b was weaker because it retained a salient
#   primary loading on the internal factor and was cleaner in promax sensitivity
#   analyses.
# - A 10-item solution that removes mdib_neg_int_remember_1b and
#   mdib_neg_ext_server_2a but retains mdib_neg_int_email_6b is therefore a
#   balanced candidate solution. A stricter 9-item solution that also removes
#   mdib_neg_int_email_6b gives the cleanest loading pattern, but leaves only
#   three internal items and does not clearly improve model fit relative to the
#   10-item solution.

# ---------------------------------------------------------------------------- #
# Check R version and load packages ----
# ---------------------------------------------------------------------------- #

source("./code/1a_define_functions.R")

groundhog_day <- version_control()

pkgs <- c("psych", "lavaan")
groundhog.library(pkgs, groundhog_day)

set.seed(1234)

# ---------------------------------------------------------------------------- #
# Define helper functions ----
# ---------------------------------------------------------------------------- #

# Create directory if it does not already exist.
make_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

# Short labels used only in output folder/file names. These keep generated paths
# short enough for collaborators to clone the repository across operating systems.
short_item_name <- function(item) {
  item_map <- c(
    mdib_neg_int_remember_1b = "1b",
    mdib_neg_ext_server_2a   = "2a",
    mdib_neg_int_email_6b    = "6b"
  )

  if (item %in% names(item_map)) {
    unname(item_map[[item]])
  } else {
    item
  }
}

# Export lavaan EFA summaries, detailed output, loadings, and the fitted object.
# The loadings file is printed output rather than a strictly rectangular CSV,
# because lavaan::efa stores loadings in an object that is most readable in print form.
export_efa_res <- function(fit, path, filename_stem) {
  make_dir(path)

  sink(file.path(path, paste0(filename_stem, ".txt")))
  print(summary(fit))
  sink()

  sink(file.path(path, paste0(filename_stem, "_detail.txt")))
  print(summary(fit, se = TRUE, zstat = TRUE, pvalue = TRUE))
  sink()

  sink(file.path(path, paste0(filename_stem, "_loadings.csv")))
  print(fit$loadings)
  sink()

  saveRDS(fit, file.path(path, paste0(filename_stem, ".rds")))
}

# Export ordinal item distributions as counts and percentages.
# These tables document response-category sparsity and provide support for
# treating the MDIB items as ordered categorical indicators.
export_item_distributions <- function(df, path, filename_stem) {
  make_dir(path)

  dist_list <- lapply(names(df), function(item) {
    tab <- table(df[[item]], useNA = "ifany")
    data.frame(
      item = item,
      response = names(tab),
      n = as.integer(tab),
      percent = round(100 * as.integer(tab) / nrow(df), 2),
      stringsAsFactors = FALSE
    )
  })

  dist_df <- do.call(rbind, dist_list)
  write.csv(
    dist_df,
    file.path(path, paste0(filename_stem, "_item_distributions.csv")),
    row.names = FALSE
  )

  invisible(dist_df)
}

# Plot item histograms for visual inspection.
# The x-axis is fixed to the 0 to 4 MDIB response scale.
plot_item_hists <- function(df, path, filename_stem, n_per_page = 6) {
  make_dir(path)

  pdf(file.path(path, paste0(filename_stem, "_hists.pdf")), height = 6, width = 6)

  for (start_col in seq(1, ncol(df), by = n_per_page)) {
    end_col <- min(start_col + n_per_page - 1, ncol(df))
    cols <- start_col:end_col

    par(mfrow = c(3, 2))
    for (i in cols) {
      hist(
        df[[i]],
        main = names(df)[i],
        xlab = "",
        breaks = seq(-0.5, 4.5, by = 1),
        xaxt = "n"
      )
      axis(1, at = 0:4)
    }
  }

  dev.off()
}

# Convert MDIB items to ordered factors for categorical analyses.
# Unobserved response categories are not artificially added. Because Script 2
# restricts the EFA sample to complete baseline MDIB data, this function should
# receive item data with no missing responses.

make_ordered_mdib <- function(df) {
  stopifnot(sum(is.na(df)) == 0)
  stopifnot(sum(df == 99, na.rm = TRUE) == 0)
  
  df_ord <- as.data.frame(lapply(df, function(x) {
    factor(x, levels = sort(unique(x)), ordered = TRUE)
  }))
  
  stopifnot(ncol(df_ord) == ncol(df))
  stopifnot(all(names(df_ord) == names(df)))
  
  df_ord
}

# Run parallel analysis based on principal components and polychoric correlations.
# Important: psych::fa.parallel prints ncomp, but for this project we do not base
# the decision on the printed ncomp value. Instead, following prior guidance for
# related analyses, we explicitly count how many observed principal-component
# eigenvalues exceed the mean simulated or resampled eigenvalues.
#
# Expected warnings:
# - "The items do not have an equal number of response alternatives" can occur
#   because severely skewed items have sparse or unobserved response categories.
# - "Matrix was not positive definite, smoothing was done" can occur for the
#   36-item polychoric matrix, likely because the number of items is large relative
#   to the usable sample size and several response categories are sparse.
run_pa_poly <- function(df, path, filename_stem, fm = "minres", n_iter = 100) {
  make_dir(path)

  pdf(file.path(path, paste0(filename_stem, "_pa_poly_", fm, ".pdf")),
      height = 6, width = 6)

  result <- psych::fa.parallel(
    df,
    fa = "pc",
    fm = fm,
    n.iter = n_iter,
    correct = 0,
    cor = "poly"
  )

  dev.off()

  n_random_mean <- sum(result$pc.values > result$pc.sim)
  n_resample_mean <- sum(result$pc.values > result$pc.simr)

  summary_df <- data.frame(
    filename_stem = filename_stem,
    fm = fm,
    n_iter = n_iter,
    threshold = c("random_data_mean", "resampled_data_mean"),
    n_components = c(n_random_mean, n_resample_mean),
    stringsAsFactors = FALSE
  )

  write.csv(
    summary_df,
    file.path(path, paste0(filename_stem, "_pa_poly_", fm, "_summary.csv")),
    row.names = FALSE
  )

  saveRDS(
    result,
    file.path(path, paste0(filename_stem, "_pa_poly_", fm, ".rds"))
  )

  list(result = result, summary = summary_df)
}

# Summarize two polychoric parallel-analysis runs and define candidate factor counts.
# If the random-data and resampled-data thresholds differ, the larger number is
# treated as the upper bound. Candidate EFA solutions then follow the preregistered
# +/- 1 approach.
summarize_pa_decision <- function(pa_minres, pa_ml, path, filename_stem) {
  make_dir(path)

  pa_summary <- rbind(pa_minres$summary, pa_ml$summary)

  upper_bound <- max(pa_summary$n_components, na.rm = TRUE)
  candidate_nfactors <- seq(max(1, upper_bound - 1), upper_bound + 1)

  decision_df <- data.frame(
    filename_stem = filename_stem,
    upper_bound = upper_bound,
    candidate_nfactors = paste(candidate_nfactors, collapse = ", "),
    stringsAsFactors = FALSE
  )

  write.csv(
    pa_summary,
    file.path(path, paste0(filename_stem, "_pa_poly_combined_summary.csv")),
    row.names = FALSE
  )

  write.csv(
    decision_df,
    file.path(path, paste0(filename_stem, "_pa_decision.csv")),
    row.names = FALSE
  )

  list(
    pa_summary = pa_summary,
    decision = decision_df,
    candidate_nfactors = candidate_nfactors
  )
}

# Run WLSMV EFAs with primary and sensitivity rotations.
# Oblimin is the primary rotation. Geomin and promax are used to evaluate whether
# the substantive loading pattern is robust to rotation choice.
run_wlsmv_efas <- function(df_ord, nfactors, path, filename_stem) {
  make_dir(path)

  rotations <- c("oblimin", "geomin", "promax")
  fits <- list()

  for (rotation in rotations) {
    set.seed(1234)

    fit <- lavaan::efa(
      data = df_ord,
      nfactors = nfactors,
      rotation = rotation,
      estimator = "WLSMV",
      ordered = names(df_ord),
      check.vcov = FALSE
    )

    fits[[rotation]] <- fit

    export_efa_res(
      fit = fit,
      path = path,
      filename_stem = paste0(filename_stem, "_", rotation, "_wlsmv")
    )
  }

  fits
}

# Run one item-removal sequence.
# Each step removes one additional item, reruns polychoric parallel analysis, and
# reruns WLSMV EFAs. These results are used to compare strategically chosen
# removal paths, not to exhaustively search every possible item subset.
run_removal_sequence <- function(df, sequence_name, removal_order, base_path) {
  sequence_path <- file.path(base_path, sequence_name)
  make_dir(sequence_path)

  current_df <- df
  sequence_log <- data.frame(
    step = integer(),
    removed_item = character(),
    retained_n_items = integer(),
    retained_items = character(),
    stringsAsFactors = FALSE
  )

  for (step in seq_along(removal_order)) {
    item_to_remove <- removal_order[step]

    if (!item_to_remove %in% names(current_df)) {
      stop(paste0("Item not found in current data: ", item_to_remove))
    }

    current_df <- current_df[, names(current_df) != item_to_remove, drop = FALSE]

    step_name <- paste0("s", sprintf("%02d", step), "_minus_", short_item_name(item_to_remove))
    step_path <- file.path(sequence_path, step_name)

    # Parallel analysis is rerun after each item removal because removing one
    # item can change the observed and simulated eigenvalue comparison.
    pa_minres <- run_pa_poly(
      current_df,
      path = file.path(step_path, "pa"),
      filename_stem = step_name,
      fm = "minres",
      n_iter = 100
    )

    pa_ml <- run_pa_poly(
      current_df,
      path = file.path(step_path, "pa"),
      filename_stem = step_name,
      fm = "ml",
      n_iter = 100
    )

    pa_decision <- summarize_pa_decision(
      pa_minres,
      pa_ml,
      path = file.path(step_path, "pa"),
      filename_stem = step_name
    )

    df_ord <- make_ordered_mdib(current_df)

    run_wlsmv_efas(
      df_ord = df_ord,
      nfactors = pa_decision$candidate_nfactors,
      path = file.path(step_path, "efa"),
      filename_stem = step_name
    )

    sequence_log <- rbind(
      sequence_log,
      data.frame(
        step = step,
        removed_item = item_to_remove,
        retained_n_items = ncol(current_df),
        retained_items = paste(names(current_df), collapse = ", "),
        stringsAsFactors = FALSE
      )
    )

    save(
      current_df,
      file = file.path(step_path, paste0(step_name, "_numeric.RData"))
    )

    save(
      df_ord,
      file = file.path(step_path, paste0(step_name, "_ordered.RData"))
    )
  }

  write.csv(
    sequence_log,
    file.path(sequence_path, paste0(sequence_name, "_sequence_log.csv")),
    row.names = FALSE
  )

  invisible(sequence_log)
}

# ---------------------------------------------------------------------------- #
# Import PD data ----
# ---------------------------------------------------------------------------- #

# The data preparation scripts should export a PD-labeled object. If this object is
# not available, revise the data preparation scripts first rather than continuing
# to use HD-labeled object names for the PD analysis.
load("./data/further_clean/mdib_pd_dat.RData")
load("./data/helper/mdib_dat_items.RData")
load("./data/helper/mdib_item_map.RData")

stopifnot(exists("mdib_pd_dat"))
stopifnot(exists("mdib_dat_items"))
stopifnot(exists("mdib_item_map"))

# ---------------------------------------------------------------------------- #
# Prepare baseline MDIB item data ----
# ---------------------------------------------------------------------------- #

mdib_items <- c(mdib_dat_items$mdib_ben, mdib_dat_items$mdib_neg)

mdib_bl <- mdib_pd_dat[
  mdib_pd_dat$redcap_event_name == "baseline_arm_1",
  mdib_items
]

# Order columns by meaning and then domain, using the predefined item map.
mdib_item_map <- mdib_item_map[order(mdib_item_map$meaning, mdib_item_map$domain), ]

mdib_bl <- mdib_bl[
  match(mdib_item_map$items_rename, names(mdib_bl))
]

stopifnot(ncol(mdib_bl) == 36)
stopifnot(all(names(mdib_bl) == mdib_item_map$items_rename))

# Confirm that the imported baseline MDIB item data are complete, as required
# by the preregistered EFA analysis sample. Therefore, no item-level imputation
# is performed in this script.

# Confirm the expected EFA analysis sample size for the current cleaned PD data
# export. Update expected_efa_n only if the cleaned data export or preregistered
# exclusion rule changes after team review.
expected_efa_n <- 82

stopifnot(nrow(mdib_bl) == length(unique(mdib_pd_dat$record_id)))
stopifnot(sum(is.na(mdib_bl)) == 0)
stopifnot(sum(mdib_bl == 99, na.rm = TRUE) == 0)
stopifnot(nrow(mdib_bl) == expected_efa_n)


# ---------------------------------------------------------------------------- #
# Define output paths ----
# ---------------------------------------------------------------------------- #

efa_path <- "./results/efa_pd/"
make_dir(efa_path)

# ---------------------------------------------------------------------------- #
# Step 1: Inspect item distributions ----
# ---------------------------------------------------------------------------- #

dist_path <- file.path(efa_path, "dist")
make_dir(dist_path)

export_item_distributions(
  df = mdib_bl,
  path = dist_path,
  filename_stem = "all36"
)

plot_item_hists(
  df = mdib_bl,
  path = dist_path,
  filename_stem = "all36"
)

# Result note:
#
# Visual inspection of the updated item distributions continued to support the
# ordered-categorical EFA workflow. The negative items showed pronounced floor
# effects, with an average of approximately 78% of responses in categories 0 or 1
# and only approximately 10% of responses in categories 3 or 4. In contrast, the
# benign items showed more responses in the middle-to-upper categories, with an
# average of approximately 46% of responses in categories 3 or 4. These
# distributional patterns support using polychoric correlations in the parallel
# analyses and WLSMV estimation in the EFAs.

# ---------------------------------------------------------------------------- #
# Step 2: Parallel analysis for all 36 MDIB items ----
# ---------------------------------------------------------------------------- #

all_items_path <- file.path(efa_path, "all36")
all_pa_path <- file.path(all_items_path, "pa")

pa_all_minres <- run_pa_poly(
  df = mdib_bl,
  path = all_pa_path,
  filename_stem = "all36",
  fm = "minres",
  n_iter = 100
)

pa_all_ml <- run_pa_poly(
  df = mdib_bl,
  path = all_pa_path,
  filename_stem = "all36",
  fm = "ml",
  n_iter = 100
)

pa_all_decision <- summarize_pa_decision(
  pa_minres = pa_all_minres,
  pa_ml = pa_all_ml,
  path = all_pa_path,
  filename_stem = "all36"
)

# Result note:
# In the current PD run, the minres and ML polychoric parallel analyses both
# supported an upper bound of three components for the 36 MDIB items. Following
# the preregistered +/- 1 approach, the candidate WLSMV EFA solutions are
# therefore 2, 3, and 4 factors.
#
# Warning note:
# The 36-item polychoric parallel analyses may produce warnings that items do not
# have the same number of observed response alternatives and that the polychoric
# correlation matrix is not positive definite. These warnings are consistent with
# sparse response categories and the complexity of estimating a 36-item
# polychoric matrix in the current usable sample.

# ---------------------------------------------------------------------------- #
# Step 3: WLSMV EFAs for all 36 MDIB items ----
# ---------------------------------------------------------------------------- #

mdib_bl_ord <- make_ordered_mdib(mdib_bl)

fits_all_36 <- run_wlsmv_efas(
  df_ord = mdib_bl_ord,
  nfactors = pa_all_decision$candidate_nfactors,
  path = file.path(all_items_path, "efa"),
  filename_stem = "all36"
)

# Result note:
# In the current PD run, model fit improved as the number of factors increased
# across the 2-, 3-, and 4-factor 36-item solutions. However, the loading patterns
# did not yield a clear, stable, theoretically interpretable structure across all
# 36 items. The 2-factor solution primarily separated benign and negative items
# rather than internal and external threat items. The 3- and 4-factor solutions
# were more complex and less stable across rotations. For this reason, the
# 36-item analyses are treated as preliminary full-item diagnostics, and the
# primary item-retention work proceeds with the 12 theorized negative items.

# ---------------------------------------------------------------------------- #
# Step 4: Restrict to the 12 theorized negative bias items ----
# ---------------------------------------------------------------------------- #

mdib_bl_neg_12 <- mdib_bl[, grepl("^mdib_neg", names(mdib_bl)), drop = FALSE]

stopifnot(ncol(mdib_bl_neg_12) == 12)

# The 12 negative items are a subset of the complete baseline MDIB item data and
# should therefore also have no missing item-level responses.
stopifnot(sum(is.na(mdib_bl_neg_12)) == 0)
stopifnot(sum(mdib_bl_neg_12 == 99, na.rm = TRUE) == 0)
stopifnot(nrow(mdib_bl_neg_12) == expected_efa_n)


neg_12_path <- file.path(efa_path, "neg12")
neg_12_dist_path <- file.path(neg_12_path, "dist")

export_item_distributions(
  df = mdib_bl_neg_12,
  path = neg_12_dist_path,
  filename_stem = "neg12"
)

plot_item_hists(
  df = mdib_bl_neg_12,
  path = neg_12_dist_path,
  filename_stem = "neg12"
)

# Result note:
# The 12 negative items show pronounced floor effects, with sparse endorsement of
# high response categories. This distributional pattern supports retaining the
# categorical-item workflow for the negative-item EFAs.

# ---------------------------------------------------------------------------- #
# Step 5: Parallel analysis for the 12 negative bias items ----
# ---------------------------------------------------------------------------- #

neg_12_pa_path <- file.path(neg_12_path, "pa")

pa_neg_12_minres <- run_pa_poly(
  df = mdib_bl_neg_12,
  path = neg_12_pa_path,
  filename_stem = "neg12",
  fm = "minres",
  n_iter = 100
)

pa_neg_12_ml <- run_pa_poly(
  df = mdib_bl_neg_12,
  path = neg_12_pa_path,
  filename_stem = "neg12",
  fm = "ml",
  n_iter = 100
)

pa_neg_12_decision <- summarize_pa_decision(
  pa_minres = pa_neg_12_minres,
  pa_ml = pa_neg_12_ml,
  path = neg_12_pa_path,
  filename_stem = "neg12"
)

# Result note:
# In the current PD run, the minres and ML polychoric parallel analyses both
# supported a one-component upper bound for the 12 negative items. Following the
# preregistered +/- 1 approach, the candidate WLSMV EFA solutions are therefore
# 1 and 2 factors.
#
# Warning note:
# The 12-negative-item polychoric parallel analyses may produce warnings that
# items do not have the same number of observed response alternatives. This is
# expected given the floor effects and sparse upper response categories.

# ---------------------------------------------------------------------------- #
# Step 6: WLSMV EFAs for the 12 negative bias items ----
# ---------------------------------------------------------------------------- #

mdib_bl_neg_12_ord <- make_ordered_mdib(mdib_bl_neg_12)

fits_neg_12 <- run_wlsmv_efas(
  df_ord = mdib_bl_neg_12_ord,
  nfactors = pa_neg_12_decision$candidate_nfactors,
  path = file.path(neg_12_path, "efa"),
  filename_stem = "neg12"
)

# Save the complete baseline 12-negative-item data used in the EFA. These
# objects should reflect the preregistered analysis sample restriction.

stopifnot(nrow(mdib_bl_neg_12) == expected_efa_n)
stopifnot(sum(is.na(mdib_bl_neg_12)) == 0)
stopifnot(sum(mdib_bl_neg_12 == 99, na.rm = TRUE) == 0)


save(
  mdib_bl_neg_12,
  file = "./data/further_clean/mdib_bl_neg_12_pd.RData"
)

save(
  mdib_bl_neg_12_ord,
  file = "./data/further_clean/mdib_bl_neg_12_ord_pd.RData"
)

# Result note:
# The 1-factor solution supported a broad negative bias factor, but absolute
# fit was poor and mdib_neg_int_remember_1b had very low communality and did
# not load saliently.
#
# The 2-factor solution improved relative fit and was more consistent with the
# theorized internal/external distinction. Most external threat items loaded on
# one factor and most internal threat items loaded on the other. However,
# absolute fit remained poor, the factor correlation was high, and key item-level
# concerns remained. mdib_neg_int_remember_1b did not load saliently,
# mdib_neg_ext_server_2a loaded primarily with the internal items despite being
# theorized as an external item, and mdib_neg_int_email_6b showed a mild
# cross-loading in the geomin solution, although its primary loading was on the
# internal factor and the oblimin/promax solutions were cleaner.
# ---------------------------------------------------------------------------- #
# Step 7: Strategically chosen item-removal sequences for negative items ----
# ---------------------------------------------------------------------------- #

# These sequences are not exhaustive. They are chosen to reflect decision rules
# that can be described in the paper:
#
# Sequence A: Start with the item that did not load saliently on either factor,
# then remove the internal item with a mild cross-loading and the external item
# with a theory-inconsistent loading.
#
# Sequence B: Start with the external item that loaded with the internal factor,
# then remove the internal item with a mild cross-loading, and then remove the
# nonsalient item if still needed.
#
# Sequence C: Start with the internal item with a mild cross-loading, then remove
# the external item that loaded with the internal factor, and then remove the
# nonsalient item if still needed.
#
# Result-informed rationale:
# Across these sequences, mdib_neg_int_remember_1b and mdib_neg_ext_server_2a
# were the clearest and most consistent problematic items. mdib_neg_int_email_6b
# was less clearly problematic because it retained a salient primary loading on
# the internal factor and was cleaner under promax rotation.

removal_sequences <- list(
  ns_int_ext = c(
    "mdib_neg_int_remember_1b",
    "mdib_neg_int_email_6b",
    "mdib_neg_ext_server_2a"
  ),
  ext_first = c(
    "mdib_neg_ext_server_2a",
    "mdib_neg_int_email_6b",
    "mdib_neg_int_remember_1b"
  ),
  int_first = c(
    "mdib_neg_int_email_6b",
    "mdib_neg_ext_server_2a",
    "mdib_neg_int_remember_1b"
  )
)

removal_path <- file.path(efa_path, "rm_seq")

for (sequence_name in names(removal_sequences)) {
  run_removal_sequence(
    df = mdib_bl_neg_12,
    sequence_name = sequence_name,
    removal_order = removal_sequences[[sequence_name]],
    base_path = removal_path
  )
}

# Result note:
# All three sequences eventually support the same strict 9-item solution after
# removing mdib_neg_int_remember_1b, mdib_neg_ext_server_2a, and
# mdib_neg_int_email_6b. This 9-item solution gives the cleanest loading pattern,
# with the retained external items loading on one factor and the retained internal
# items loading on the other. However, this strict solution leaves only three
# internal items and does not clearly improve model fit relative to the more
# balanced 10-item solution that retains mdib_neg_int_email_6b.

# ---------------------------------------------------------------------------- #
# Step 8: Focused sequence for comparing the 10-item and 9-item candidates ----
# ---------------------------------------------------------------------------- #

# This focused sequence starts with the two most consistently problematic items:
# - mdib_neg_int_remember_1b: nonsalient loading and very low communality.
# - mdib_neg_ext_server_2a: external item that repeatedly loaded with internal items.
#
# Step 2 of this sequence gives the balanced 10-item candidate:
# remove mdib_neg_int_remember_1b and mdib_neg_ext_server_2a, retain
# mdib_neg_int_email_6b.
#
# Step 3 gives the stricter 9-item candidate:
# also remove mdib_neg_int_email_6b.

removal_sequences <- list(
  focused = c(
    "mdib_neg_int_remember_1b",
    "mdib_neg_ext_server_2a",
    "mdib_neg_int_email_6b"
  )
)

removal_path <- file.path(efa_path, "rm_focus")

for (sequence_name in names(removal_sequences)) {
  run_removal_sequence(
    df = mdib_bl_neg_12,
    sequence_name = sequence_name,
    removal_order = removal_sequences[[sequence_name]],
    base_path = removal_path
  )
}

# Final interpretation note:
# Based on the current PD EFA results, the 10-item solution that removes
# mdib_neg_int_remember_1b and mdib_neg_ext_server_2a but retains
# mdib_neg_int_email_6b appears to be the
# most balanced reduced-item candidate. This solution removes the two clearest
# problematic items, preserves four internal items, and yields a
# clear internal/external two-factor pattern. The 9-item solution that also removes
# mdib_neg_int_email_6b can be retained as a stricter sensitivity or alternative
# solution because it gives the cleanest loading pattern but leaves only three
# internal items.
#
# These comments document the decision logic for discussion with the team. The
# final retained item set should be decided after reviewing the EFA notes, factor
# loading tables, model fit, item content, and theoretical coverage.

# ---------------------------------------------------------------------------- #
# End of script ----
# ---------------------------------------------------------------------------- #
