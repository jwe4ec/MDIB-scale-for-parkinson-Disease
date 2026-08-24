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
# This script excludes prior analyses that were run but were vestiges of the HD script:
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
#   component upper bound. Following the preregistered +/- 1 approach and theory
#   positing 2 factors, 1- and 2-factor WLSMV EFA solutions were examined.
# - The 12-item 1-factor solution supported a broad negative bias factor but had
#   poor model fit and a very weak item, mdib_neg_int_remember_1b.
# - The 12-item 2-factor solution was more consistent with the theorized
#   internal/external distinction, but three item-level concerns remained:
#     1. mdib_neg_int_remember_1b did not load saliently on either factor.
#     2. mdib_neg_ext_server_2a loaded with the internal rather than external items.
#     3. mdib_neg_int_email_6b showed a mild cross-loading in some rotations.
# - Across the item-removal sequences, mdib_neg_int_remember_1b and
#   mdib_neg_ext_server_2a were the clearest problematic items. The evidence for
#   removing mdib_neg_int_email_6b was weaker because it retained a salient primary 
#   loading on the internal factor and was cleaner in promax sensitivity analyses.
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
source("./code/4a_define_efa_functions.R")

groundhog_day <- version_control()

pkgs <- c("psych", "lavaan")
groundhog.library(pkgs, groundhog_day)

set.seed(1234)

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
stopifnot(all(range(mdib_bl, na.rm = TRUE) == c(0, 4)))
stopifnot(sum(is.na(mdib_bl)) == 0)
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

# Result notes:
# - Items with empty categories:
#   - mdib_ben_int_remember_1c: no responses of 0
#   - mdib_neg_ext_server_2a: no responses of 4
# - Visual inspection of the updated item distributions continued to support the
#   ordered-categorical EFA workflow. The negative items showed pronounced floor
#   effects, with an average of approximately 78% of responses in categories 0 or 1
#   and only approximately 10% of responses in categories 3 or 4. By contrast, the
#   benign items showed more responses in the middle-to-upper categories, with an
#   average of approximately 46% of responses in categories 3 or 4. These
#   distributional patterns support using polychoric correlations in the parallel
#   analyses and WLSMV estimation in the EFAs.

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

export_item_distributions(
  df = mdib_bl_neg_12,
  path = dist_path,
  filename_stem = "neg12"
)

plot_item_hists(
  df = mdib_bl_neg_12,
  path = dist_path,
  filename_stem = "neg12"
)

# Result note:
# As noted above, the 12 negative items show pronounced floor effects, with 
# sparse endorsement of high response categories. This distributional pattern 
# supports retaining the categorical-item workflow for the negative-item EFAs.

# ---------------------------------------------------------------------------- #
# Step 5: Parallel analysis for the 12 negative bias items ----
# ---------------------------------------------------------------------------- #

neg_12_path <- file.path(efa_path, "neg12")
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
# preregistered +/- 1 approach (and theory), the candidate WLSMV EFA solutions 
# are therefore 1 and 2 factors.
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

# Save the complete baseline 12-negative-item data used in the EFA

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
# theorized as an external item, and mdib_neg_int_email_6b showed a moderate
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