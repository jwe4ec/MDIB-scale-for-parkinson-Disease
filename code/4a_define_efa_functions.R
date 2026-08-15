# ---------------------------------------------------------------------------- #
# Define Helper Functions for Exploratory Factor Analyses
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
# (2_clean_data_compute_item_missingness.R) restricts the EFA sample to complete 
# baseline MDIB data, this function should receive item data with no missing responses.
make_ordered_mdib <- function(df) {
  stopifnot(all(range(df, na.rm = TRUE) == c(0, 4)))
  stopifnot(sum(is.na(df)) == 0)
  
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
# the HD analyses, we explicitly count how many observed principal-component
# eigenvalues exceed the mean of the eigenvalues across many datasets (based on
# either simulated random data or resampling from our sample's data).
#
# Expected warnings:
# - "In polychoric(x, correct = correct): The items do not have an equal number 
#   of response alternatives, global set to FALSE."
#   - Can occur because severely skewed items have sparse or unobserved response categories.
# - "In cor.smooth(mat): Matrix was not positive definite, smoothing was done"
#   - Can occur for the 36-item polychoric matrix, likely because the number of 
#     items is large relative to the usable sample size and several response 
#     categories are sparse.
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