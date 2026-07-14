# ---------------------------------------------------------------------------- #
# Store working directory, load helper functions, and set package-version date ----
# ---------------------------------------------------------------------------- #

# Store the project root directory. This object can be used later if the script
# needs to return to the original working directory after writing outputs.

wd_dir <- getwd()

# Load custom helper functions used across the analysis scripts, including
# version_control().

source("./code/1a_define_functions.R")

# Check the R version used for the reproducibility snapshot, load groundhog, and
# store the package-version date used when loading analysis packages.

groundhog_day <- version_control()

# ---------------------------------------------------------------------------- #
# Import PD Aim 1 data ----
# ---------------------------------------------------------------------------- #

# Import the de-identified PD Aim 1 REDCap data.

mdib_pd_dat <- read.csv("./data/bot_cleaned/final PD Aim 1 data_deid_2022-12-08_OSF.csv")

# ---------------------------------------------------------------------------- #
# Define REDCap event names ----
# ---------------------------------------------------------------------------- #

# Define baseline and follow-up event names used throughout the data cleaning and
# scoring scripts.

bl   <- "baseline_arm_1"
fu   <- "followup_arm_1"
both <- c(bl, fu)

# ---------------------------------------------------------------------------- #
# Check for blank rows ----
# ---------------------------------------------------------------------------- #

# Confirm that the current PD data export does not contain blank rows with missing
# record_id values.

stopifnot(sum(is.na(mdib_pd_dat$record_id)) == 0)

# ---------------------------------------------------------------------------- #
# Remove participants who consented but never started the baseline survey ----
# ---------------------------------------------------------------------------- #

# Identify survey item columns, excluding participant ID and REDCap event name.

target_cols <- setdiff(names(mdib_pd_dat), c("record_id", "redcap_event_name"))

# Identify baseline rows for participants who consented but never started the
# survey. These rows contain only NA, 0, or "" across all survey item columns.

row_never_started <- vector(length = nrow(mdib_pd_dat))

for (i in 1:nrow(mdib_pd_dat)) {
  row_never_started[[i]] <- mdib_pd_dat$redcap_event_name[[i]] == bl &
    all(apply(mdib_pd_dat[i, target_cols], 2, function(x) {
      is.na(x) | x %in% c(0, "")
    }))
}

# Confirm that 7 participants meet this criterion.

stopifnot(sum(row_never_started) == 7)
stopifnot(length(unique(mdib_pd_dat$record_id[row_never_started])) == 7)

# Remove participants who consented but never started the baseline survey.

mdib_pd_dat <- mdib_pd_dat[!row_never_started, ]

# Confirm that 88 unique participants remain.

stopifnot(length(unique(mdib_pd_dat$record_id)) == 88)

# ---------------------------------------------------------------------------- #
# Identify items for relevant scales ----
# ---------------------------------------------------------------------------- #

# Identify item variables for scales used in the PD Aim 1 analyses, including
# negative and benign interpretation bias items from the MDIB and BBSIQ, anxiety
# sensitivity items from the ASI, fear of negative evaluation items from the
# BFNE-II, anxiety symptoms from Neuro-QoL Anxiety, social avoidance and distress
# items from the full SADS at baseline and reduced SADS at follow-up, and alcohol
# use items from the AUDIT-C. The identified item lists are checked below against
# the expected variable names in the current PD data export.

nms <- names(mdib_pd_dat)

mdib_neg_items <- nms[startsWith(nms, "md_bbsiq_") & endsWith(nms, "neg")]
mdib_ben_items <- nms[startsWith(nms, "md_bbsiq_") & endsWith(nms, "benign")]

# BBSIQ items are identified separately from MDIB items by excluding variables
# that begin with the MDIB-specific "md_bbsiq" prefix.

bbsiq_neg_items_mdib <- nms[startsWith(nms, "bbsiq_") & endsWith(nms, "neg")]
bbsiq_ben_items_mdib <- nms[startsWith(nms, "bbsiq_") & endsWith(nms, "benign")]

asi_items <- nms[startsWith(nms, "asi_")]

bfne2_items <- nms[startsWith(nms, "bfne_")]

neuroqol_anx_items <- nms[startsWith(nms, "neuroqol_") & !endsWith(nms, "_complete")]

sads_items     <- nms[startsWith(nms, "sad_") & !endsWith(nms, "_v2")]
sads_red_items <- nms[startsWith(nms, "sad_") & endsWith(nms, "_v2")]

auditc_items <- nms[startsWith(nms, "alcohol_audit_c")]

# Confirm that the identified item variables match the expected variable names in
# the current PD data export.

stopifnot(identical(
  mdib_neg_items,
  c("md_bbsiq_1b_neg", "md_bbsiq_2a_neg", "md_bbsiq_3c_neg",
    "md_bbsiq_4c_neg", "md_bbsiq_5a_neg", "md_bbsiq_6b_neg",
    "md_bbsiq_7a_neg", "md_bbsiq_8b_neg", "md_bbsiq_9c_neg",
    "md_bbsiq_10a_neg", "md_bbsiq_11c_neg", "md_bbsiq_12b_neg")
))

stopifnot(identical(
  mdib_ben_items,
  c("md_bbsiq_1a_benign", "md_bbsiq_1c_benign", "md_bbsiq_2b_benign",
    "md_bbsiq_2c_benign", "md_bbsiq_3a_benign", "md_bbsiq_3b_benign",
    "md_bbsiq_4a_benign", "md_bbsiq_4b_benign", "md_bbsiq_5b_benign",
    "md_bbsiq_5c_benign", "md_bbsiq_6a_benign", "md_bbsiq_6c_benign",
    "md_bbsiq_7b_benign", "md_bbsiq_7c_benign", "md_bbsiq_8a_benign",
    "md_bbsiq_8c_benign", "md_bbsiq_9a_benign", "md_bbsiq_9b_benign",
    "md_bbsiq_10b_benign", "md_bbsiq_10c_benign", "md_bbsiq_11a_benign",
    "md_bbsiq_11b_benign", "md_bbsiq_12a_benign", "md_bbsiq_12c_benign")
))

stopifnot(identical(
  bbsiq_neg_items_mdib,
  c("bbsiq_1c_neg", "bbsiq_2b_neg", "bbsiq_3c_neg",
    "bbsiq_4c_neg", "bbsiq_5a_neg", "bbsiq_6a_neg",
    "bbsiq_7b_neg", "bbsiq_8c_neg", "bbsiq_9b_neg",
    "bbsiq_10b_neg", "bbsiq_11b_neg", "bbsiq_12a_neg",
    "bbsiq_13c_neg", "bbsiq_14c_neg")
))

stopifnot(identical(
  bbsiq_ben_items_mdib,
  c("bbsiq_1a_benign", "bbsiq_1b_benign", "bbsiq_2a_benign",
    "bbsiq_2c_benign", "bbsiq_3a_benign", "bbsiq_3b_benign",
    "bbsiq_4a_benign", "bbsiq_4b_benign", "bbsiq_5b_benign",
    "bbsiq_5c_benign", "bbsiq_6b_benign", "bbsiq_6c_benign",
    "bbsiq_7a_benign", "bbsiq_7c_benign", "bbsiq_8a_benign",
    "bbsiq_8b_benign", "bbsiq_9a_benign", "bbsiq_9c_benign",
    "bbsiq_10a_benign", "bbsiq_10c_benign", "bbsiq_11a_benign",
    "bbsiq_11c_benign", "bbsiq_12b_benign", "bbsiq_12c_benign",
    "bbsiq_13a_benign", "bbsiq_13b_benign", "bbsiq_14a_benign",
    "bbsiq_14b_benign")
))

stopifnot(identical(
  asi_items,
  c("asi_1", "asi_2", "asi_3", "asi_4", "asi_5", "asi_6", "asi_7",
    "asi_8", "asi_9", "asi_10", "asi_11", "asi_12", "asi_13", "asi_14",
    "asi_15", "asi_16")
))

stopifnot(identical(
  bfne2_items,
  c("bfne_1", "bfne_2", "bfne_3", "bfne_4", "bfne_5", "bfne_6",
    "bfne_7", "bfne_8", "bfne_9", "bfne_10", "bfne_11", "bfne_12")
))

stopifnot(identical(
  neuroqol_anx_items,
  c("neuroqol_edanx53", "neuroqol_edanx46", "neuroqol_edanx48",
    "neuroqol_edanx41", "neuroqol_edanx54", "neuroqol_edanx55",
    "neuroqol_edanx18", "neuroqol_nqanx07")
))

stopifnot(identical(
  sads_items,
  c("sad_1", "sad_2", "sad_3", "sad_4", "sad_5", "sad_6", "sad_7",
    "sad_8", "sad_9", "sad_10", "sad_11", "sad_12", "sad_13", "sad_14",
    "sad_15", "sad_16", "sad_17", "sad_18", "sad_19", "sad_20", "sad_21",
    "sad_22", "sad_23", "sad_24", "sad_25", "sad_26", "sad_27", "sad_28")
))

stopifnot(identical(
  sads_red_items,
  c("sad_20_v2", "sad_27_v2", "sad_13_v2", "sad_12_v2",
    "sad_24_v2", "sad_15_v2", "sad_4_v2", "sad_16_v2")
))

stopifnot(identical(
  auditc_items,
  c("alcohol_audit_c_1", "alcohol_audit_c_2", "alcohol_audit_c_3")
))

# Confirm the expected number of items for each scale.

stopifnot(length(mdib_neg_items) == 12)
stopifnot(length(mdib_ben_items) == 24)
stopifnot(length(bbsiq_neg_items_mdib) == 14)
stopifnot(length(bbsiq_ben_items_mdib) == 28)
stopifnot(length(asi_items) == 16)
stopifnot(length(bfne2_items) == 12)
stopifnot(length(neuroqol_anx_items) == 8)
stopifnot(length(sads_items) == 28)
stopifnot(length(sads_red_items) == 8)
stopifnot(length(auditc_items) == 3)

# ---------------------------------------------------------------------------- #
# Rename MDIB items ----
# ---------------------------------------------------------------------------- #

# Rename MDIB items so that variable names indicate interpretation valence
# (benign vs. negative), threat domain (internal vs. external), scenario, and
# original item number. Store a crosswalk between the original REDCap variable
# names and the renamed variables for reproducibility.

mdib_items <- c(mdib_ben_items, mdib_neg_items)

mdib_items_rename <-
  c("mdib_ben_int_remember_1a", "mdib_ben_int_remember_1c", "mdib_ben_ext_server_2b",     "mdib_ben_ext_server_2c",
    "mdib_ben_ext_reminder_3a", "mdib_ben_ext_reminder_3b", "mdib_ben_int_cleaning_4a",   "mdib_ben_int_cleaning_4b",
    "mdib_ben_ext_neighbor_5b", "mdib_ben_ext_neighbor_5c", "mdib_ben_int_email_6a",      "mdib_ben_int_email_6c",
    "mdib_ben_ext_exercise_7b", "mdib_ben_ext_exercise_7c", "mdib_ben_int_medication_8a", "mdib_ben_int_medication_8c",
    "mdib_ben_ext_walk_9a",     "mdib_ben_ext_walk_9b",     "mdib_ben_ext_job_10b",       "mdib_ben_ext_job_10c",
    "mdib_ben_ext_stumble_11a", "mdib_ben_ext_stumble_11b", "mdib_ben_int_cough_12a",     "mdib_ben_int_cough_12c",
    "mdib_neg_int_remember_1b", "mdib_neg_ext_server_2a",   "mdib_neg_ext_reminder_3c",   "mdib_neg_int_cleaning_4c",
    "mdib_neg_ext_neighbor_5a", "mdib_neg_int_email_6b",    "mdib_neg_ext_exercise_7a",   "mdib_neg_int_medication_8b",
    "mdib_neg_ext_walk_9c",     "mdib_neg_ext_job_10a",     "mdib_neg_ext_stumble_11c",   "mdib_neg_int_cough_12b")

# Extract item numbers from the renamed MDIB variables and pad one-digit item
# numbers with a leading zero so that items sort correctly.

item_number <- sub(".*_(\\d+[a-z])$", "\\1", mdib_items_rename) 
item_number <- sub("^([0-9])([a-z])$", "0\\1\\2", item_number)  

# Define item-level metadata used to create the MDIB item crosswalk.

meaning <-
  c("ben", "ben", "ben", "ben",
    "ben", "ben", "ben", "ben",
    "ben", "ben", "ben", "ben",
    "ben", "ben", "ben", "ben",
    "ben", "ben", "ben", "ben",
    "ben", "ben", "ben", "ben",
    "neg", "neg", "neg", "neg",
    "neg", "neg", "neg", "neg",
    "neg", "neg", "neg", "neg")

domain <-
  c("int", "int", "ext", "ext",
    "ext", "ext", "int", "int",
    "ext", "ext", "int", "int",
    "ext", "ext", "int", "int",
    "ext", "ext", "ext", "ext",
    "ext", "ext", "int", "int",
    "int", "ext", "ext", "int",
    "ext", "int", "ext", "int",
    "ext", "ext", "ext", "int")

# Create item-level crosswalk with the original item name, renamed item name,
# item number, interpretation valence, and threat domain.

mdib_item_map <- data.frame(items        = mdib_items,
                            items_rename = mdib_items_rename,
                            item_number  = item_number,
                            meaning      = meaning,
                            domain       = domain)

# Apply the renamed MDIB item variables to the PD data.

names(mdib_pd_dat)[names(mdib_pd_dat) %in% mdib_items] <- 
  mdib_items_rename[match(names(mdib_pd_dat)[names(mdib_pd_dat) %in% mdib_items], mdib_items)]

# Update the MDIB item lists to use the renamed variables.

mdib_neg_items <- mdib_item_map$items_rename[mdib_item_map$meaning == "neg"]
mdib_ben_items <- mdib_item_map$items_rename[mdib_item_map$meaning == "ben"]

# Confirm the expected number of renamed negative and benign MDIB items.

stopifnot(length(mdib_neg_items) == 12)
stopifnot(length(mdib_ben_items) == 24)

# ---------------------------------------------------------------------------- #
# Define scale items ----
# ---------------------------------------------------------------------------- #

# Define theory-based MDIB negative-bias item sets by threat domain. Internal
# threat items reflect catastrophizing about disease progression, and external
# threat items reflect negative social evaluation.

mdib_neg_int_items <- mdib_item_map$items_rename[
  mdib_item_map$meaning == "neg" &
    mdib_item_map$domain == "int"
]

mdib_neg_ext_items <- mdib_item_map$items_rename[
  mdib_item_map$meaning == "neg" &
    mdib_item_map$domain == "ext"
]

# Confirm that the theory-based MDIB negative-bias item sets match the expected
# internal and external threat items.

stopifnot(identical(
  mdib_neg_int_items,
  c("mdib_neg_int_remember_1b",
    "mdib_neg_int_cleaning_4c",
    "mdib_neg_int_email_6b",
    "mdib_neg_int_medication_8b",
    "mdib_neg_int_cough_12b")
))

stopifnot(identical(
  mdib_neg_ext_items,
  c("mdib_neg_ext_server_2a",
    "mdib_neg_ext_reminder_3c",
    "mdib_neg_ext_neighbor_5a",
    "mdib_neg_ext_exercise_7a",
    "mdib_neg_ext_walk_9c",
    "mdib_neg_ext_job_10a",
    "mdib_neg_ext_stumble_11c")
))

stopifnot(length(mdib_neg_int_items) == 5)
stopifnot(length(mdib_neg_ext_items) == 7)

# Define reduced MDIB negative-bias item sets retained in the current PD EFA
# solution. These item sets should be updated only if the team revises the final
# retained item set after reviewing the alternative item-exclusion sequences.

mdib_neg_9_int_items <- c(
  "mdib_neg_int_cleaning_4c",
  "mdib_neg_int_medication_8b",
  "mdib_neg_int_cough_12b"
)

mdib_neg_9_ext_items <- c(
  "mdib_neg_ext_reminder_3c",
  "mdib_neg_ext_neighbor_5a",
  "mdib_neg_ext_exercise_7a",
  "mdib_neg_ext_walk_9c",
  "mdib_neg_ext_job_10a",
  "mdib_neg_ext_stumble_11c"
)

stopifnot(length(mdib_neg_9_int_items) == 3)
stopifnot(length(mdib_neg_9_ext_items) == 6)

# Define theory-based BBSIQ negative-bias item sets by threat domain. The suffix
# "_mdib" is retained because these are the BBSIQ item names as they appear in
# the MDIB/PD REDCap dataset, not because they refer to the HD sample or HD
# analysis.

bbsiq_neg_int_items_mdib <- c(
  "bbsiq_2b_neg",
  "bbsiq_3c_neg",
  "bbsiq_5a_neg",
  "bbsiq_8c_neg",
  "bbsiq_11b_neg",
  "bbsiq_12a_neg",
  "bbsiq_14c_neg"
)

bbsiq_neg_ext_items_mdib <- c(
  "bbsiq_1c_neg",
  "bbsiq_4c_neg",
  "bbsiq_6a_neg",
  "bbsiq_7b_neg",
  "bbsiq_9b_neg",
  "bbsiq_10b_neg",
  "bbsiq_13c_neg"
)

stopifnot(length(bbsiq_neg_int_items_mdib) == 7)
stopifnot(length(bbsiq_neg_ext_items_mdib) == 7)

# Define reduced ASI subscale item sets for physical, cognitive, and social
# concerns. These reduced item sets are based on prior three-factor solutions and
# item-retention decisions described in the ASI scoring documentation.

asi_red_phy_items <- c("asi_11", "asi_9")
asi_red_cog_items <- c("asi_12", "asi_2")
asi_red_soc_items <- c("asi_1", "asi_13")

stopifnot(length(asi_red_phy_items) == 2)
stopifnot(length(asi_red_cog_items) == 2)
stopifnot(length(asi_red_soc_items) == 2)

# Define the preferred 8-item BFNE-II scale by removing the four reverse-scored
# items excluded by Carleton et al. (2007).

bfne2_8_items <- bfne2_items[
  !(bfne2_items %in% c("bfne_2", "bfne_4", "bfne_7", "bfne_11"))
]

stopifnot(identical(
  bfne2_8_items,
  c("bfne_1", "bfne_3", "bfne_5", "bfne_6",
    "bfne_8", "bfne_9", "bfne_10", "bfne_12")
))

stopifnot(length(bfne2_8_items) == 8)

# Store all item sets in a single helper list used for scoring, missing-data
# summaries, and later analysis scripts.

mdib_dat_items <- list(
  mdib_neg       = mdib_neg_items,
  mdib_neg_int   = mdib_neg_int_items,
  mdib_neg_ext   = mdib_neg_ext_items,
  mdib_neg_9_int = mdib_neg_9_int_items,
  mdib_neg_9_ext = mdib_neg_9_ext_items,
  mdib_ben       = mdib_ben_items,
  bbsiq_neg      = bbsiq_neg_items_mdib,
  bbsiq_neg_int  = bbsiq_neg_int_items_mdib,
  bbsiq_neg_ext  = bbsiq_neg_ext_items_mdib,
  bbsiq_ben      = bbsiq_ben_items_mdib,
  asi            = asi_items,
  asi_red_phy    = asi_red_phy_items,
  asi_red_cog    = asi_red_cog_items,
  asi_red_soc    = asi_red_soc_items,
  bfne2          = bfne2_items,
  bfne2_8        = bfne2_8_items,
  neuroqol_anx   = neuroqol_anx_items,
  sads           = sads_items,
  sads_red       = sads_red_items,
  auditc         = auditc_items
)

# ---------------------------------------------------------------------------- #
# Inspect item-level NA values before recoding planned missing-value codes ----
# ---------------------------------------------------------------------------- #

# Inspect item-level NA values in the current PD data export before recoding
# planned missing-value codes and computing scale scores. Missingness is
# summarized by scale within the time points where each scale is used.

# Define helper function for counting item-level NA values within specified
# REDCap event names.

count_item_na <- function(dat, items, time_points) {
  dat_sub <- dat[dat$redcap_event_name %in% time_points, items]
  sum(is.na(dat_sub))
}


# Summarize item-level NA values by scale.

na_diagnostic_tbl <- data.frame(
  scale = c("mdib_neg", "mdib_ben", "bbsiq_neg", "bbsiq_ben",
            "neuroqol_anx", "asi", "bfne2", "sads"),
  time_points = c("baseline_followup", "baseline_followup",
                  "baseline_followup", "baseline_followup",
                  "baseline_followup", "baseline", "baseline", "baseline"),
  n_na = c(
    count_item_na(mdib_pd_dat, mdib_dat_items$mdib_neg, both),
    count_item_na(mdib_pd_dat, mdib_dat_items$mdib_ben, both),
    count_item_na(mdib_pd_dat, mdib_dat_items$bbsiq_neg, both),
    count_item_na(mdib_pd_dat, mdib_dat_items$bbsiq_ben, both),
    count_item_na(mdib_pd_dat, mdib_dat_items$neuroqol_anx, both),
    count_item_na(mdib_pd_dat, mdib_dat_items$asi, bl),
    count_item_na(mdib_pd_dat, mdib_dat_items$bfne2, bl),
    count_item_na(mdib_pd_dat, mdib_dat_items$sads, bl)
  )
)

# Confirm the expected item-level NA counts in the current PD data export.

stopifnot(
  all(na_diagnostic_tbl$n_na == c(81, 162, 114, 228, 64, 64, 48, 112))
)

#na_diagnostic_tbl


# ---------------------------------------------------------------------------- #
# Inspect planned AUDIT-C skip patterns ----
# ---------------------------------------------------------------------------- #

# AUDIT-C has planned NA values at baseline because REDCap skipped items based on
# alcohol-use screening responses. These skipped items are recoded to 0 later in
# the script when they indicate no lifetime alcohol use or no current alcohol use.

auditc_bl <- mdib_pd_dat[
  mdib_pd_dat$redcap_event_name == bl,
  c("record_id", "alcohol_ever", mdib_dat_items$auditc)
]

n_obs_na_auditc <- sum(is.na(auditc_bl[, mdib_dat_items$auditc]))
stopifnot(n_obs_na_auditc == 85)

# REDCap skipped all 3 AUDIT-C items when alcohol_ever was 0, indicating no
# lifetime alcohol use.

n_na_auditc_alcohol_never <- sum(is.na(
  auditc_bl[auditc_bl$alcohol_ever == 0, mdib_dat_items$auditc]
))
stopifnot(n_na_auditc_alcohol_never == 51)

# REDCap skipped AUDIT-C Items 2 and 3 when Item 1 was 0, indicating no current
# alcohol use.

n_na_auditc_item1_never <- sum(is.na(
  auditc_bl[
    !is.na(auditc_bl$alcohol_audit_c_1) &
      auditc_bl$alcohol_audit_c_1 == 0,
    c("alcohol_audit_c_2", "alcohol_audit_c_3")
  ]
))
stopifnot(n_na_auditc_item1_never == 34)

# ---------------------------------------------------------------------------- #
# Inspect planned reduced SADS non-administration at follow-up ----
# ---------------------------------------------------------------------------- #

# Reduced SADS has planned NA values at follow-up because the reduced SADS was
# not administered to 5 participants.

sads_red_fu <- mdib_pd_dat[
  mdib_pd_dat$redcap_event_name == fu,
  c("record_id", mdib_dat_items$sads_red)
]

n_obs_na_sads_red <- sum(is.na(sads_red_fu[, mdib_dat_items$sads_red]))
stopifnot(n_obs_na_sads_red == 40)

rows_all_items_na_sads_red <- rowSums(
  !is.na(sads_red_fu[, mdib_dat_items$sads_red])
) == 0

n_rows_all_items_na_sads_red <- sum(rows_all_items_na_sads_red)
stopifnot(n_rows_all_items_na_sads_red == 5)

stopifnot(
  n_rows_all_items_na_sads_red * length(mdib_dat_items$sads_red) ==
    n_obs_na_sads_red
)


# ---------------------------------------------------------------------------- #
# Recode planned AUDIT-C skip values ----
# ---------------------------------------------------------------------------- #

# Recode planned AUDIT-C skip values as 0 at baseline. REDCap skipped AUDIT-C
# items when participants reported no lifetime alcohol use or no current alcohol
# use. These skipped items indicate non-use rather than item nonresponse.

# If AUDIT-C Item 1 was 0 ("never"), Items 2 and 3 were skipped because they were
# not applicable. Recode Items 2 and 3 as 0 for these participants.

auditc_item1_never_bl <- mdib_pd_dat$redcap_event_name == bl &
  !is.na(mdib_pd_dat$alcohol_audit_c_1) &
  mdib_pd_dat$alcohol_audit_c_1 == 0

mdib_pd_dat[
  auditc_item1_never_bl,
  c("alcohol_audit_c_2", "alcohol_audit_c_3")
] <- 0

# If alcohol_ever was 0, all 3 AUDIT-C items were skipped because the participant
# reported no lifetime alcohol use. Recode all 3 AUDIT-C items as 0 for these
# participants.

alcohol_never_bl <- mdib_pd_dat$redcap_event_name == bl &
  !is.na(mdib_pd_dat$alcohol_ever) &
  mdib_pd_dat$alcohol_ever == 0

mdib_pd_dat[
  alcohol_never_bl,
  mdib_dat_items$auditc
] <- 0

# Confirm that the planned AUDIT-C skip values have been recoded.

stopifnot(
  sum(is.na(mdib_pd_dat[
    alcohol_never_bl,
    mdib_dat_items$auditc
  ])) == 0
)

stopifnot(
  sum(is.na(mdib_pd_dat[
    auditc_item1_never_bl,
    c("alcohol_audit_c_2", "alcohol_audit_c_3")
  ])) == 0
)


# ---------------------------------------------------------------------------- #
# Identify participants with incomplete MDIB data at baseline ----
# ---------------------------------------------------------------------------- #

# The preregistered analysis sample is restricted to participants with complete
# item-level MDIB data at baseline. Incomplete baseline MDIB data are defined as
# having at least one baseline MDIB item that is either truly missing (NA) or
# coded as 99 ("prefer not to answer").
#
# Earlier versions of this script only flagged baseline MDIB items coded as 99.
# That was too narrow because some participants had NA values on baseline MDIB
# items and were therefore retained in the exported analysis object. The code
# below flags both sources of incomplete baseline item-level MDIB data.

mdib_items <- c(mdib_dat_items$mdib_neg, mdib_dat_items$mdib_ben)

mdib_bl_item_dat <- mdib_pd_dat[
  mdib_pd_dat$redcap_event_name == bl,
  c("record_id", mdib_items)
]

mdib_bl_item_mat <- mdib_bl_item_dat[, mdib_items]

rows_incomplete_mdib_bl <- rowSums(
  is.na(mdib_bl_item_mat) | mdib_bl_item_mat == 99,
  na.rm = TRUE
) > 0

incompl_mdib_bl_data_ids <- mdib_bl_item_dat$record_id[
  rows_incomplete_mdib_bl
]

# Create a diagnostic table showing why each participant is excluded from the
# baseline MDIB complete-case analysis sample.

incompl_mdib_bl_tbl <- data.frame(
  record_id = incompl_mdib_bl_data_ids,
  n_mdib_na_bl = rowSums(is.na(mdib_bl_item_mat))[rows_incomplete_mdib_bl],
  n_mdib_pna_bl = rowSums(mdib_bl_item_mat == 99, na.rm = TRUE)[rows_incomplete_mdib_bl],
  stringsAsFactors = FALSE
)

incompl_mdib_bl_tbl$n_mdib_incomplete_bl <-
  incompl_mdib_bl_tbl$n_mdib_na_bl + incompl_mdib_bl_tbl$n_mdib_pna_bl

# Confirm that the identified IDs are unique and that every retained participant
# has complete baseline MDIB item-level data.

stopifnot(length(incompl_mdib_bl_data_ids) == length(unique(incompl_mdib_bl_data_ids)))
stopifnot(all(incompl_mdib_bl_tbl$n_mdib_incomplete_bl > 0))

# In the current PD data export, the complete-case baseline MDIB restriction
# removes 6 participants from the 88 participants who started the baseline survey,
# leaving 82 participants in the analysis sample. If these checks fail after a
# future data export, inspect incompl_mdib_bl_tbl and update the expected counts.

stopifnot(length(incompl_mdib_bl_data_ids) == 6)
stopifnot(88 - length(incompl_mdib_bl_data_ids) == 82)

# Export the diagnostic table for reproducibility.

missing_rates_path <- "./results/missingness/"
dir.create(missing_rates_path, recursive = TRUE, showWarnings = FALSE)

write.csv(
  incompl_mdib_bl_tbl,
  file.path(missing_rates_path, "incomplete_mdib_bl_exclusion_tbl.csv"),
  row.names = FALSE
)


# ---------------------------------------------------------------------------- #
# Compute scale-level missingness due to all items coded as "prefer not to answer" ----
# ---------------------------------------------------------------------------- #

# Before recoding "prefer not to answer" values as NA, summarize the number of
# scale scores that will be missing because all items in a given scale are coded
# as 99. Participants with incomplete baseline MDIB data are excluded from this
# summary because they are removed from the analysis sample later in this script.

temp_dat <- mdib_pd_dat[
  !(mdib_pd_dat$record_id %in% incompl_mdib_bl_data_ids),
]

# Define helper function to count the number of scale scores across specified
# REDCap events for which all items are coded as 99.

compute_all_item_missingness <- function(dat, scale, items, time_points) {
  dat <- dat[dat$redcap_event_name %in% time_points, ]
  
  pna <- 99
  
  rows_all_items_pna <- rowSums(dat[, items] == pna, na.rm = TRUE) == length(items)
  n_rows_all_items_pna <- sum(rows_all_items_pna)
  
  cat(scale, ": ", n_rows_all_items_pna, "\n", sep = "")
  print(table(dat[rows_all_items_pna, "redcap_event_name"]))
  cat("\n", "-----", "\n\n")
}

# Write all-item PNA missingness summaries to a text file.

missing_rates_path <- "./results/missingness/"
dir.create(missing_rates_path, recursive = TRUE, showWarnings = FALSE)

sink(file = paste0(missing_rates_path, "all_item_missingness.txt"))

cat("Number of Scale Scores Missing Due to 'Prefer Not to Answer' for All Items:", "\n\n")

compute_all_item_missingness(temp_dat, "mdib_neg_9_int_m", mdib_dat_items$mdib_neg_9_int, both)
compute_all_item_missingness(temp_dat, "mdib_neg_9_ext_m", mdib_dat_items$mdib_neg_9_ext, both)
compute_all_item_missingness(temp_dat, "bbsiq_neg_int_m",  mdib_dat_items$bbsiq_neg_int,  both)
compute_all_item_missingness(temp_dat, "bbsiq_neg_ext_m",  mdib_dat_items$bbsiq_neg_ext,  both)
compute_all_item_missingness(temp_dat, "asi_m",            mdib_dat_items$asi,            bl)
compute_all_item_missingness(temp_dat, "asi_red_phy_m",    mdib_dat_items$asi_red_phy,    bl)
compute_all_item_missingness(temp_dat, "asi_red_cog_m",    mdib_dat_items$asi_red_cog,    bl)
compute_all_item_missingness(temp_dat, "asi_red_soc_m",    mdib_dat_items$asi_red_soc,    bl)
compute_all_item_missingness(temp_dat, "bfne2_8_m",        mdib_dat_items$bfne2_8,        bl)
compute_all_item_missingness(temp_dat, "neuroqol_anx_m",   mdib_dat_items$neuroqol_anx,   both)
compute_all_item_missingness(temp_dat, "sads_m",           mdib_dat_items$sads,           bl)
compute_all_item_missingness(temp_dat, "sads_red_m",       mdib_dat_items$sads_red,       fu)
compute_all_item_missingness(temp_dat, "auditc_m",         mdib_dat_items$auditc,         bl)

sink()

# ---------------------------------------------------------------------------- #
# Recode "prefer not to answer" values ----
# ---------------------------------------------------------------------------- #

# Recode "prefer not to answer" values, coded as 99, as NA for all item-level
# variables used in scale scoring.

target_items <- c(
  mdib_dat_items$mdib_neg,
  mdib_dat_items$mdib_ben,
  mdib_dat_items$bbsiq_neg,
  mdib_dat_items$bbsiq_ben,
  mdib_dat_items$asi,
  mdib_dat_items$bfne2,
  mdib_dat_items$neuroqol_anx,
  mdib_dat_items$sads,
  mdib_dat_items$sads_red,
  mdib_dat_items$auditc
)

target_dat <- mdib_pd_dat[, target_items]
target_dat[target_dat == 99] <- NA
mdib_pd_dat[, target_items] <- target_dat

# ---------------------------------------------------------------------------- #
# Recode BFNE-II items ----
# ---------------------------------------------------------------------------- #

# Recode BFNE-II items from the REDCap 1-5 response scale to the 0-4 scale used
# by Carleton et al. (2007). The response options shown in REDCap did not include
# numeric labels, so subtracting 1 aligns the stored values with the scoring scale.

stopifnot(all(range(mdib_pd_dat[, mdib_dat_items$bfne2], na.rm = TRUE) == c(1, 5)))

mdib_pd_dat[, mdib_dat_items$bfne2] <- mdib_pd_dat[, mdib_dat_items$bfne2] - 1

stopifnot(all(range(mdib_pd_dat[, mdib_dat_items$bfne2], na.rm = TRUE) == c(0, 4)))

# ---------------------------------------------------------------------------- #
# Compute scale scores ----
# ---------------------------------------------------------------------------- #

# Compute scale scores as the mean of available items. For MDIB negative bias,
# use the reduced item sets from the current PD EFA solution. These item sets
# should be updated if the team revises the final retained MDIB items.

compute_row_mean <- function(dat, items) {
  score <- rowMeans(dat[, items], na.rm = TRUE)
  score[is.nan(score)] <- NA
  score
}

mdib_pd_dat$mdib_neg_9_int_m <- compute_row_mean(mdib_pd_dat, mdib_dat_items$mdib_neg_9_int)
mdib_pd_dat$mdib_neg_9_ext_m <- compute_row_mean(mdib_pd_dat, mdib_dat_items$mdib_neg_9_ext)
mdib_pd_dat$bbsiq_neg_int_m  <- compute_row_mean(mdib_pd_dat, mdib_dat_items$bbsiq_neg_int)
mdib_pd_dat$bbsiq_neg_ext_m  <- compute_row_mean(mdib_pd_dat, mdib_dat_items$bbsiq_neg_ext)
mdib_pd_dat$asi_m            <- compute_row_mean(mdib_pd_dat, mdib_dat_items$asi)
mdib_pd_dat$asi_red_phy_m    <- compute_row_mean(mdib_pd_dat, mdib_dat_items$asi_red_phy)
mdib_pd_dat$asi_red_cog_m    <- compute_row_mean(mdib_pd_dat, mdib_dat_items$asi_red_cog)
mdib_pd_dat$asi_red_soc_m    <- compute_row_mean(mdib_pd_dat, mdib_dat_items$asi_red_soc)
mdib_pd_dat$bfne2_8_m        <- compute_row_mean(mdib_pd_dat, mdib_dat_items$bfne2_8)
mdib_pd_dat$neuroqol_anx_m   <- compute_row_mean(mdib_pd_dat, mdib_dat_items$neuroqol_anx)
mdib_pd_dat$sads_m           <- compute_row_mean(mdib_pd_dat, mdib_dat_items$sads)
mdib_pd_dat$sads_red_m       <- compute_row_mean(mdib_pd_dat, mdib_dat_items$sads_red)
mdib_pd_dat$auditc_m         <- compute_row_mean(mdib_pd_dat, mdib_dat_items$auditc)

# ---------------------------------------------------------------------------- #
# Create table of item-level missingness for MDIB at baseline ----
# ---------------------------------------------------------------------------- #

# Summarize item-level missingness for all MDIB items at baseline before removing
# participants with incomplete baseline MDIB data.

mdib_bl <- mdib_pd_dat[
  mdib_pd_dat$redcap_event_name == bl,
  c(mdib_dat_items$mdib_neg, mdib_dat_items$mdib_ben)
]

# Order MDIB items by original item number.

mdib_item_map <- mdib_item_map[order(mdib_item_map$item_number), ]

mdib_bl <- mdib_bl[
  match(mdib_item_map$items_rename, names(mdib_bl))
]

# Confirm the baseline MDIB missingness table is based on 88 participants before
# removing participants with incomplete baseline MDIB data.

n_mdib_bl <- nrow(mdib_bl)
stopifnot(n_mdib_bl == 88)

item_n_missing <- colSums(is.na(mdib_bl))
item_perc_missing <- format(
  round((item_n_missing / n_mdib_bl) * 100, 1),
  nsmall = 1,
  trim = TRUE
)
item_n_perc_missing <- paste0(item_n_missing, " (", item_perc_missing, ")")

mdib_bl_item_missing_tbl <- data.frame(
  domain         = mdib_item_map$domain,
  item           = names(mdib_bl),
  meaning        = mdib_item_map$meaning,
  n_missing      = item_n_missing,
  n_perc_missing = item_n_perc_missing
)

# Retain only MDIB items with at least one missing baseline response.

mdib_bl_item_missing_tbl <- mdib_bl_item_missing_tbl[
  mdib_bl_item_missing_tbl$n_missing > 0,
]

row.names(mdib_bl_item_missing_tbl) <- 1:nrow(mdib_bl_item_missing_tbl)

# Export item-level baseline MDIB missingness table.

write.csv(
  mdib_bl_item_missing_tbl,
  paste0(missing_rates_path, "mdib_bl_item_missing_tbl.csv"),
  row.names = FALSE
)

# ---------------------------------------------------------------------------- #
# Remove participants with incomplete MDIB data at baseline ----
# ---------------------------------------------------------------------------- #

# Remove participants with incomplete baseline MDIB data, defined above as
# having at least one baseline MDIB item that is either NA or coded as 99
# ("prefer not to answer"). This implements the preregistered complete-case
# baseline MDIB item-level restriction and leaves an analysis sample of 82
# participants in the current PD data export.

mdib_pd_dat <- mdib_pd_dat[
  !(mdib_pd_dat$record_id %in% incompl_mdib_bl_data_ids),
]

stopifnot(length(unique(mdib_pd_dat$record_id)) == 82)

mdib_bl_after_exclusion <- mdib_pd_dat[
  mdib_pd_dat$redcap_event_name == bl,
  mdib_items
]

stopifnot(sum(is.na(mdib_bl_after_exclusion)) == 0)
stopifnot(sum(mdib_bl_after_exclusion == 99, na.rm = TRUE) == 0)

# ---------------------------------------------------------------------------- #
# Compute rates of item-level missingness in computed scale scores ----
# ---------------------------------------------------------------------------- #

# For each computed scale score, calculate the percentage of observed scale
# scores that were computed with at least one missing item but not all items
# missing.

compute_some_item_missingness <- function(dat, scale, items, time_points) {
  dat <- dat[dat$redcap_event_name %in% time_points, ]
  
  denom <- sum(!is.na(dat[, scale]))
  
  rows_at_least_one_item_na <- rowSums(is.na(dat[, items])) > 0
  rows_all_items_na <- rowSums(!is.na(dat[, items])) == 0
  
  numer <- nrow(dat[rows_at_least_one_item_na & !rows_all_items_na, ])
  
  prop <- numer / denom
  percent <- prop * 100
  
  cat(scale, ": ", percent, "%", "\n", sep = "")
}

sink(file = paste0(missing_rates_path, "some_item_missingness.txt"))

cat("Percentages of Scale Scores Computed With At Least One Item Missing:", "\n\n")

compute_some_item_missingness(mdib_pd_dat, "mdib_neg_9_int_m", mdib_dat_items$mdib_neg_9_int, both)
compute_some_item_missingness(mdib_pd_dat, "mdib_neg_9_ext_m", mdib_dat_items$mdib_neg_9_ext, both)
compute_some_item_missingness(mdib_pd_dat, "bbsiq_neg_int_m",  mdib_dat_items$bbsiq_neg_int,  both)
compute_some_item_missingness(mdib_pd_dat, "bbsiq_neg_ext_m",  mdib_dat_items$bbsiq_neg_ext,  both)
compute_some_item_missingness(mdib_pd_dat, "asi_m",            mdib_dat_items$asi,            bl)
compute_some_item_missingness(mdib_pd_dat, "asi_red_phy_m",    mdib_dat_items$asi_red_phy,    bl)
compute_some_item_missingness(mdib_pd_dat, "asi_red_cog_m",    mdib_dat_items$asi_red_cog,    bl)
compute_some_item_missingness(mdib_pd_dat, "asi_red_soc_m",    mdib_dat_items$asi_red_soc,    bl)
compute_some_item_missingness(mdib_pd_dat, "bfne2_8_m",        mdib_dat_items$bfne2_8,        bl)
compute_some_item_missingness(mdib_pd_dat, "neuroqol_anx_m",   mdib_dat_items$neuroqol_anx,   both)
compute_some_item_missingness(mdib_pd_dat, "sads_m",           mdib_dat_items$sads,           bl)
compute_some_item_missingness(mdib_pd_dat, "sads_red_m",       mdib_dat_items$sads_red,       fu)
compute_some_item_missingness(mdib_pd_dat, "auditc_m",         mdib_dat_items$auditc,         bl)

sink()

# ---------------------------------------------------------------------------- #
# Export cleaned data and helper objects ----
# ---------------------------------------------------------------------------- #

# Create output folders if they do not already exist.

dir.create("./data/further_clean", recursive = TRUE, showWarnings = FALSE)
dir.create("./data/helper", recursive = TRUE, showWarnings = FALSE)

# Export cleaned PD data and helper objects for later analysis scripts.

save(mdib_pd_dat, file = "./data/further_clean/mdib_pd_dat.RData")
save(mdib_dat_items, file = "./data/helper/mdib_dat_items.RData")
save(mdib_item_map,  file = "./data/helper/mdib_item_map.RData")








