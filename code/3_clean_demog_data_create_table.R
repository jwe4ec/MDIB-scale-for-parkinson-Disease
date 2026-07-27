# ---------------------------------------------------------------------------- #
# Clean PD demographic data and create demographics table
# ---------------------------------------------------------------------------- #

# This script imports the cleaned PD Aim 1 analysis data, extracts baseline
# demographic variables, cleans and labels demographic variables, creates
# demographics tables, and exports CSV and Word versions of the tables.

# Before running this script, restart R and set the working directory to the
# project root folder.

# ---------------------------------------------------------------------------- #
# Load helper functions and set package-version date ----
# ---------------------------------------------------------------------------- #

source("./code/1a_define_functions.R")

groundhog_day <- version_control()

pkgs <- c("flextable", "officer", "ftExtra")
groundhog.library(pkgs, groundhog_day)

source("./code/1b_set_flextable_defaults.R")
source("./code/1c_set_officer_properties.R")

# ---------------------------------------------------------------------------- #
# Import cleaned PD analysis data ----
# ---------------------------------------------------------------------------- #

load("./data/further_clean/mdib_pd_dat.RData")

bl <- "baseline_arm_1"

stopifnot(exists("mdib_pd_dat"))
stopifnot(length(unique(mdib_pd_dat$record_id)) == 82)

# ---------------------------------------------------------------------------- #
# Extract baseline demographic data ----
# ---------------------------------------------------------------------------- #

# Extract demographic variables from the baseline PD analysis sample.

index_cols <- c("record_id", "redcap_event_name")

dem_cols <- c(
  "age",
  paste0("race___", c(1:5, 9, 99)),
  "race_other",
  "ethnicity",
  "country",
  "relationship_status",
  "live_alone",
  "gender",
  "sex",
  "employment_status",
  "education",
  "age_diagnosis_pd",
  "survey_help",
  paste0("study_awareness___", c(0:5, 8))
)

dem_dat <- mdib_pd_dat[
  mdib_pd_dat$redcap_event_name == bl,
  c(index_cols, dem_cols)
]

# Confirm that the demographics table is based on the 82 participants in the
# baseline PD analysis sample.

stopifnot(nrow(dem_dat) == 82)
stopifnot(length(unique(dem_dat$record_id)) == 82)

# Store a copy of the raw baseline demographic data before cleaning.

dem_raw <- dem_dat

# ---------------------------------------------------------------------------- #
# Clean age and age at PD diagnosis ----
# ---------------------------------------------------------------------------- #

# Confirm missingness and observed range for age in the baseline PD analysis
# sample.

stopifnot(sum(is.na(dem_dat$age)) == 2)
stopifnot(all(range(dem_dat$age, na.rm = TRUE) == c(40, 84)))

# Confirm missingness and observed range for age at PD diagnosis in the baseline
# PD analysis sample.

stopifnot(sum(is.na(dem_dat$age_diagnosis_pd)) == 2)
stopifnot(all(range(dem_dat$age_diagnosis_pd, na.rm = TRUE) == c(35, 79)))

# ---------------------------------------------------------------------------- #
# Clean race ----
# ---------------------------------------------------------------------------- #

# One participant entered "Latino" in the race-other text field. Because Latino is
# an ethnicity response rather than a race response, clear the race-other text
# field and remove the corresponding "Other" race checkbox. This participant
# also selected other race categories and is categorized as "More than one race"
# in the collapsed race variable below.

dem_dat$race_other <- trimws(dem_dat$race_other)

race_other_latino <- grepl(
  "^latino$",
  dem_dat$race_other,
  ignore.case = TRUE
)

stopifnot(sum(race_other_latino, na.rm = TRUE) == 1)

dem_dat$race_other[race_other_latino] <- ""
dem_dat$race___9[race_other_latino] <- 0

# Confirm the number of participants who selected "Other" race and "Prefer not to
# answer" for race after correcting the Latino race-other entry.

stopifnot(sum(dem_dat$race___9, na.rm = TRUE) == 0)
stopifnot(sum(dem_dat$race___99, na.rm = TRUE) == 1)

# Confirm that the race-other text field is blank after cleaning.

stopifnot(all(dem_dat$race_other == "" | is.na(dem_dat$race_other)))

# Create a single race variable, collapsing participants who selected more than
# one non-PNA race category into "More than one race."

dem_dat$race_coll <- NA

race_cols_not_pna <- paste0("race___", c(1:5, 9))

for (i in 1:nrow(dem_dat)) {
  n_race_selected <- rowSums(dem_dat[i, race_cols_not_pna], na.rm = TRUE)
  
  if (n_race_selected <= 1) {
    dem_dat$race_coll[[i]][dem_dat$race___1[[i]]  == 1] <- "American Indian or Alaska Native"
    dem_dat$race_coll[[i]][dem_dat$race___2[[i]]  == 1] <- "Asian"
    dem_dat$race_coll[[i]][dem_dat$race___3[[i]]  == 1] <- "Black or African American"
    dem_dat$race_coll[[i]][dem_dat$race___4[[i]]  == 1] <- "Native Hawaiian or Other Pacific Islander"
    dem_dat$race_coll[[i]][dem_dat$race___5[[i]]  == 1] <- "White"
    dem_dat$race_coll[[i]][dem_dat$race___9[[i]]  == 1] <- "Other"
    dem_dat$race_coll[[i]][dem_dat$race___99[[i]] == 1] <- "Prefer not to answer"
  } else if (n_race_selected > 1) {
    dem_dat$race_coll[[i]] <- "More than one race"
  }
}

# Recode missing race values as "Missing" so that they are displayed in the
# demographics table.

dem_dat$race_coll <- as.character(dem_dat$race_coll)
dem_dat$race_coll[is.na(dem_dat$race_coll)] <- "Missing"

dem_dat$race_coll <- factor(
  dem_dat$race_coll,
  levels = c(
    "American Indian or Alaska Native",
    "Asian",
    "Black or African American",
    "Native Hawaiian or Other Pacific Islander",
    "White",
    "More than one race",
    "Other",
    "Prefer not to answer",
    "Missing"
  )
)

# Confirm the cleaned race distribution in the baseline PD analysis sample.

stopifnot(sum(dem_dat$race_coll == "Black or African American", na.rm = TRUE) == 1)
stopifnot(sum(dem_dat$race_coll == "White", na.rm = TRUE) == 75)
stopifnot(sum(dem_dat$race_coll == "More than one race", na.rm = TRUE) == 3)
stopifnot(sum(dem_dat$race_coll == "Other", na.rm = TRUE) == 0)
stopifnot(sum(dem_dat$race_coll == "Prefer not to answer", na.rm = TRUE) == 1)
stopifnot(sum(dem_dat$race_coll == "Missing", na.rm = TRUE) == 2)

table(dem_dat$race_coll, useNA = "ifany")

# ---------------------------------------------------------------------------- #
# Clean categorical demographic variables ----
# ---------------------------------------------------------------------------- #

# Recode REDCap numeric response codes into labeled factors for categorical
# demographic variables.

dem_dat$ethnicity <- factor(
  dem_dat$ethnicity,
  levels = 1:4,
  labels = c(
    "Hispanic or Latino",
    "Not Hispanic or Latino",
    "Unknown",
    "Prefer not to answer"
  )
)

dem_dat$relationship_status <- factor(
  dem_dat$relationship_status,
  levels = c(0:9, 99),
  labels = c(
    "Single",
    "Single, but casually dating",
    "Single, but currently engaged to be married",
    "Single, but currently living with someone in a marriage-like relationship",
    "Married",
    "In a domestic or civil union",
    "Separated",
    "Divorced",
    "Widow/widower",
    "Other",
    "Prefer not to answer"
  )
)

dem_dat$live_alone <- factor(
  dem_dat$live_alone,
  levels = c(0:1, 99),
  labels = c("No", "Yes", "Prefer not to answer")
)

dem_dat$gender <- factor(
  dem_dat$gender,
  levels = c(0:4, 99),
  labels = c(
    "Female",
    "Male",
    "Transgender Female",
    "Transgender Male",
    "Other",
    "Prefer not to answer"
  )
)

dem_dat$sex <- factor(
  dem_dat$sex,
  levels = c(0:3, 99),
  labels = c(
    "Female",
    "Male",
    "Intersex",
    "Unknown or other",
    "Prefer not to answer"
  )
)

dem_dat$employment_status <- factor(
  dem_dat$employment_status,
  levels = c(0:8, 99),
  labels = c(
    "Working full-time",
    "Working part-time",
    "Unemployed or laid off",
    "Looking for work",
    "Homemaker/keeping house or raising children full-time",
    "Retired",
    "Student",
    "Other",
    "Unknown",
    "Prefer not to answer"
  )
)

dem_dat$education <- factor(
  dem_dat$education,
  levels = c(1:8, 99),
  labels = c(
    "Elementary School",
    "Junior High",
    "High School",
    "Some College",
    "Associate's Degree",
    "Bachelor's Degree",
    "Master's Degree",
    "Doctorate/ PhD",
    "Prefer not to answer"
  )
)

# Recode missing values in categorical demographic variables as "Missing" while
# preserving the original factor-level order.

add_missing_level <- function(x) {
  original_levels <- levels(x)
  x <- as.character(x)
  x[is.na(x)] <- "Missing"
  factor(x, levels = c(original_levels, "Missing"))
}

categorical_dem_vars <- c(
  "ethnicity",
  "relationship_status",
  "live_alone",
  "gender",
  "sex",
  "employment_status",
  "education"
)

dem_dat[categorical_dem_vars] <- lapply(
  dem_dat[categorical_dem_vars],
  add_missing_level
)

# Confirm that each categorical demographic variable has 82 observations after
# missing values are recoded as "Missing".

stopifnot(all(sapply(
  dem_dat[categorical_dem_vars],
  function(x) sum(table(x, useNA = "ifany")) == 82
)))

# Confirm that no categorical demographic variable has remaining NA values.

stopifnot(all(sapply(
  dem_dat[categorical_dem_vars],
  function(x) sum(is.na(x)) == 0
)))

# Inspect cleaned categorical demographic distributions.

lapply(
  dem_dat[categorical_dem_vars],
  table,
  useNA = "ifany"
)

# ---------------------------------------------------------------------------- #
# Clean country ----
# ---------------------------------------------------------------------------- #

# Recode variations of United States and missing country values for the extended
# demographics table.

dem_dat$country <- trimws(dem_dat$country)

united_states <- c(
  "America", "U.S.", "united states", "United States",
  "United States of America", "us", "US", "usa", "Usa", "USA"
)

dem_dat$country[dem_dat$country %in% united_states] <- "United States"
dem_dat$country[is.na(dem_dat$country) | dem_dat$country == ""] <- "Missing"

dem_dat$country <- factor(dem_dat$country)

stopifnot(sum(table(dem_dat$country, useNA = "ifany")) == 82)
stopifnot(sum(is.na(dem_dat$country)) == 0)

table(dem_dat$country, useNA = "ifany")

# ---------------------------------------------------------------------------- #
# Clean study awareness ----
# ---------------------------------------------------------------------------- #

# Create a single study-awareness variable from REDCap checkbox indicators.
# Participants who selected more than one source are categorized as
# "More than one source." Participants who did not select any source are coded as
# "Missing" for display in the demographics table.

study_awareness_cols <- paste0("study_awareness___", c(0:5, 8))

study_awareness_labels <- c(
  "Parkinson's Pals website",
  "Michael J. Fox Foundation website or email",
  "Facebook group posting",
  "International Parkinson Society website or email",
  "My PD healthcare provider",
  "A friend or family member",
  "Other"
)

n_study_awareness_selected <- rowSums(
  dem_dat[, study_awareness_cols] == 1,
  na.rm = TRUE
)

# Confirm the number of participants who selected no source or more than one
# source.

stopifnot(sum(n_study_awareness_selected == 0) == 2)
stopifnot(sum(n_study_awareness_selected > 1) == 4)
stopifnot(max(n_study_awareness_selected, na.rm = TRUE) == 2)

# Create cleaned study-awareness variable.

dem_dat$study_awareness <- "Missing"

for (i in seq_along(study_awareness_cols)) {
  dem_dat$study_awareness[
    dem_dat[[study_awareness_cols[i]]] == 1 &
      n_study_awareness_selected == 1
  ] <- study_awareness_labels[i]
}

dem_dat$study_awareness[n_study_awareness_selected > 1] <- "More than one source"

dem_dat$study_awareness <- factor(
  dem_dat$study_awareness,
  levels = c(study_awareness_labels, "More than one source", "Missing")
)

# Confirm that the cleaned study-awareness variable covers all 82 participants.

stopifnot(sum(table(dem_dat$study_awareness, useNA = "ifany")) == 82)
stopifnot(sum(is.na(dem_dat$study_awareness)) == 0)

# Inspect cleaned study-awareness distribution.

table(dem_dat$study_awareness, useNA = "ifany")

# ---------------------------------------------------------------------------- #
# Clean survey help ----
# ---------------------------------------------------------------------------- #

# Recode survey-help response codes into labeled categories. Missing values are
# coded as "Missing" for display in the demographics table.

dem_dat$survey_help <- factor(
  dem_dat$survey_help,
  levels = c(0, 1, 99),
  labels = c("No", "Yes", "Prefer not to answer")
)

dem_dat$survey_help <- as.character(dem_dat$survey_help)
dem_dat$survey_help[is.na(dem_dat$survey_help)] <- "Missing"

dem_dat$survey_help <- factor(
  dem_dat$survey_help,
  levels = c("No", "Yes", "Prefer not to answer", "Missing")
)

# Confirm that the cleaned survey-help variable covers all 82 participants.

stopifnot(sum(table(dem_dat$survey_help, useNA = "ifany")) == 82)
stopifnot(sum(is.na(dem_dat$survey_help)) == 0)

# Inspect cleaned survey-help distribution.

table(dem_dat$survey_help, useNA = "ifany")

# ---------------------------------------------------------------------------- #
# Save cleaned demographic data ----
# ---------------------------------------------------------------------------- #

dir.create("./data/further_clean", recursive = TRUE, showWarnings = FALSE)

save(dem_dat, file = "./data/further_clean/dem_dat.RData")

# ---------------------------------------------------------------------------- #
# Define function to compute demographic descriptives ----
# ---------------------------------------------------------------------------- #

compute_dem_desc <- function(df, exclude_fct_cols = character()) {
  
  # Compute sample size.
  
  n <- data.frame(label = "n", value = nrow(df))
  
  # Compute mean and standard deviation for numeric demographic variables.
  
  num_vars_all   <- c("age", "age_diagnosis_pd")
  num_vars_exist <- intersect(num_vars_all, names(df))
  
  num_labels_map <- list(
    age              = list(var = "Age",                 unit = "Years: M (SD)"),
    age_diagnosis_pd = list(var = "Age at PD Diagnosis", unit = "Years: M (SD)")
  )
  
  num_res <- data.frame()
  
  if (length(num_vars_exist) > 0) {
    for (v in num_vars_exist) {
      m   <- mean(df[[v]], na.rm = TRUE)
      sdv <- sd(df[[v]],   na.rm = TRUE)
      
      block <- rbind(
        data.frame(label = num_labels_map[[v]]$var,  value = NA_character_),
        data.frame(label = num_labels_map[[v]]$unit, value = sprintf("%.2f (%.2f)", m, sdv))
      )
      
      num_res <- rbind(num_res, block)
    }
  }
  
  # Compute counts and percentages for categorical demographic variables.
  
  fct_vars_all <- c(
    "gender",
    "sex",
    "race_coll",
    "ethnicity",
    "education",
    "employment_status",
    "relationship_status",
    "live_alone",
    "country",
    "study_awareness",
    "survey_help"
  )
  
  fct_labels_all <- c(
    "Gender: n (%)",
    "Sex Assigned at Birth: n (%)",
    "Race: n (%)",
    "Ethnicity: n (%)",
    "Education: n (%)",
    "Employment Status: n (%)",
    "Relationship Status: n (%)",
    "Living Alone: n (%)",
    "Country: n (%)",
    "Where did you hear about this survey study?: n (%)",
    "Did anyone help you complete these surveys?: n (%)"
  )
  
  names(fct_labels_all) <- fct_vars_all
  
  fct_vars_exist <- setdiff(
    intersect(fct_vars_all, names(df)),
    exclude_fct_cols
  )
  
  fct_res <- data.frame()
  
  if (length(fct_vars_exist) > 0) {
    for (v in fct_vars_exist) {
      
      tbl <- table(df[[v]])
      tbl <- tbl[tbl > 0]
      
      prop_tbl <- prop.table(tbl) * 100
      
      block <- rbind(
        data.frame(label = fct_labels_all[[v]], value = NA_character_),
        data.frame(
          label = names(tbl),
          value = paste0(
            as.numeric(tbl),
            " (",
            sprintf("%.1f", as.numeric(prop_tbl)),
            ")"
          )
        )
      )
      
      fct_res <- rbind(fct_res, block)
    }
  }
  
  pieces <- list(n, num_res, fct_res)
  pieces <- pieces[sapply(pieces, nrow) > 0]
  
  res <- do.call(rbind, pieces)
  rownames(res) <- NULL
  
  return(res)
}

# ---------------------------------------------------------------------------- #
# Create demographics tables ----
# ---------------------------------------------------------------------------- #

exclude_fct_cols <- c("country", "study_awareness", "survey_help")

dem_tbl     <- compute_dem_desc(dem_dat, exclude_fct_cols)
dem_tbl_ext <- compute_dem_desc(dem_dat, NULL)

dem_path <- "./results/demographics/"

dir.create(dem_path, recursive = TRUE, showWarnings = FALSE)

write.csv(dem_tbl,     paste0(dem_path, "dem_tbl.csv"),          row.names = FALSE)
write.csv(dem_tbl_ext, paste0(dem_path, "dem_tbl_extended.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------- #
# Format demographics table ----
# ---------------------------------------------------------------------------- #

# Format the demographics table for export to Word. Flextable defaults and Word
# section properties are sourced above.

format_dem_tbl <- function(dem_tbl, title) {
  
  # Format the label column using Markdown.
  
  dem_tbl$label_md <- dem_tbl$label
  
  rows_no_indent <- dem_tbl$label_md == "n" |
    grepl(
      "\\b(Age|Gender|Sex Assigned at Birth|Race|Ethnicity|Education|Employment Status|Relationship Status|Living Alone)\\b",
      dem_tbl$label_md
    )
  
  rows_indent <- !rows_no_indent
  
  indent_spaces <- "\\ \\ \\ \\ \\ "
  dem_tbl$label_md[rows_indent] <- paste0(indent_spaces, dem_tbl$label_md[rows_indent])
  
  dem_tbl$label_md[dem_tbl$label_md == "n"] <- "*n*"
  dem_tbl$label_md <- gsub("n \\(%\\)", "*n* \\(%\\)", dem_tbl$label_md)
  dem_tbl$label_md <- gsub("M \\(SD\\)", "*M* \\(*SD*\\)", dem_tbl$label_md)
  
  dem_tbl <- dem_tbl[c("label_md", names(dem_tbl)[names(dem_tbl) != "label_md"])]
  
  # Define columns to display.
  
  target_cols <- names(dem_tbl)[!(names(dem_tbl) %in% "label")]
  
  # Create flextable.
  
  dem_tbl_ft <- flextable(dem_tbl[, target_cols]) |>
    set_table_properties(align = "left") |>
    set_caption(
      as_paragraph(as_i(title)),
      word_stylename = "heading 1",
      fp_p = fp_par(padding.left = 0, padding.right = 0),
      align_with_table = FALSE
    ) |>
    align(align = "center", part = "header") |>
    align(align = "center", part = "body") |>
    align(j = "label_md", align = "left", part = "body") |>
    align(align = "left", part = "footer") |>
    valign(valign = "bottom", part = "header") |>
    set_header_labels(
      label_md = "Characteristic",
      value    = "Value"
    ) |>
    colformat_md(j = "label_md", part = "body") |>
    labelizor(
      part = "body",
      labels = c(
        "Doctorate/ PhD"    = "Doctorate/PhD",
        "Working full-time" = "Working full time",
        "Working part-time" = "Working part time"
      )
    ) |>
    autofit()
  
  return(dem_tbl_ft)
}

# Create formatted demographics table.

dem_tbl_ft <- format_dem_tbl(
  dem_tbl,
  title = "Sample Characteristics at Baseline"
)

# ---------------------------------------------------------------------------- #
# Write demographics table to MS Word ----
# ---------------------------------------------------------------------------- #

# Write the formatted demographics table to a Word document.

dem_tbl_orientation <- "p"
dem_tbl_number      <- 1

doc <- read_docx()
doc <- body_set_default_section(doc, psect_prop)

doc <- body_add_fpar(
  doc,
  fpar(
    ftext(
      paste0("Table ", dem_tbl_number),
      prop = text_prop_bold
    )
  )
)

doc <- body_add_par(doc, "")

doc <- body_add_flextable(doc, dem_tbl_ft, align = "left")

if (dem_tbl_orientation == "p") {
  doc <- body_end_block_section(doc, block_section(psect_prop))
} else if (dem_tbl_orientation == "l") {
  doc <- body_end_block_section(doc, block_section(lsect_prop))
}

print(doc, target = paste0(dem_path, "dem_tbl.docx"))