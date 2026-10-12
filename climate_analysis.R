################################################################################
# Intersectionality and Departmental Climate Among STEM Doctoral Students
# Analysis script (blinded replication code)
#
# What this script does
#   1. Factor analysis (EFA) of the 35 climate/experience survey items
#        - all students  (appendix/supplementary; N = 2,486)
#        - U.S. citizens / permanent residents only (N = 1,021; used for the
#          composite climate score)
#   2. Builds the composite climate score (weighted mean of the 35 items) and
#      splits it into quartiles (Q1 = lowest ... Q4 = highest climate).
#   3. "Who experiences good climate": Gender x Racial identification x
#      Institution type (PWI / HBCU)
#        - descriptive crosstab with climate quartile (Table 1)
#        - chi-square test with standardized residuals (Table 2)
#        - three-way Aligned Rank Transform (ART) ANOVA on the climate score
#          (Table 3)
#
# Inputs  (in ./data)
#   students_all_blinded.csv        all respondents, items only
#   students_citizens_blinded.csv   citizens/permanent residents, items +
#                                   gender + racial-identification variable
#
# Outputs (written to ./output; read by manuscript_tables.Rmd)
#   Table / figure files named after the manuscript element they feed
#   (Table1_..., Table2_..., Table3_..., TableA1_..., TableA2_..., Figure1_...);
#   files starting with "Supp_" or "FigS" are supplementary.
#
# How to run
#   Open the repository folder as the working directory (or an RStudio
#   project) and run the whole script. Then knit manuscript_tables.Rmd.
#
# Notes
#   - Gender and racial-identification "Prefer not to say / Other" responses
#     (and missing values) are combined into one level and RETAINED in the
#     descriptive crosstab and chi-square test. The ART ANOVA is restricted
#     to Male/Female x URM/NotURM because the combined level leaves empty
#     cells in the HBCU group (a full factorial cannot be estimated).
#   - Institution is collapsed to PWI vs. HBCU (the InstType column).
#   - Department is included in the data files for completeness but is not
#     used in these analyses.
#   - The climate score weights come from the citizens-only EFA and are used
#     unchanged for every analysis below.
################################################################################


# ------------------------------------------------------------------------------
# 0. Setup
# ------------------------------------------------------------------------------

library(psych)        # polychoric(), fa(), fa.parallel()
library(GPArotation)  # REQUIRED for oblimin rotation inside psych::fa(); without
                      # it psych returns an UNROTATED solution and all
                      # downstream weights/scores change
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)      # fct_reorder()
library(DescTools)    # CramerV()
library(ARTool)       # Aligned Rank Transform ANOVA

# Random-number seed. Parallel analysis (fa.parallel) and the Monte Carlo
# chi-square both use random numbers, so the order of those calls below
# matters for exact reproducibility of the Monte Carlo p-value.
set.seed(123)

path_all <- file.path("data", "students_all_blinded.csv")
path_cit <- file.path("data", "students_citizens_blinded.csv")

out_dir <- "output"
if (!dir.exists(out_dir)) dir.create(out_dir)

# Read a data file. Blank cells are treated as missing (NA); check.names = FALSE
# keeps the survey item text exactly as written in the file header.
read_data <- function(path) {
  read.csv(path, check.names = FALSE, stringsAsFactors = FALSE,
           na.strings = c("", "NA"))
}

# Helper: cumulative variance accounted for by a psych::fa() solution.
# (A 1-factor solution has no "Cumulative Var" row, so fall back to the
# proportion of variance.)
get_cum_var <- function(fa_obj) {
  va <- fa_obj$Vaccounted
  rn <- rownames(va)
  if ("Cumulative Var" %in% rn) {
    va["Cumulative Var", ncol(va)]
  } else if ("Proportion Var" %in% rn) {
    sum(va["Proportion Var", ])
  } else {
    NA_real_
  }
}


# ------------------------------------------------------------------------------
# 1. Factor analysis
# ------------------------------------------------------------------------------

# Columns that are NOT survey items. Everything else is one of the 35 items,
# in questionnaire order.
id_cols_all <- c("Institution", "InstType", "Department")
id_cols_cit <- c(id_cols_all, "Gender", "USURM")

# ---- 1.1 All students: 8-factor EFA (reported in Methods; supplementary) -----

dataall    <- read_data(path_all)
item_names <- setdiff(names(dataall), id_cols_all)
stopifnot(length(item_names) == 35)
dataallord <- dataall[, item_names]

# Polychoric correlations are appropriate for ordinal (Likert-type) items.
# (psych notes that items have different numbers of response options -- 5-point
#  vs. 4-point -- and so does not assume a common 'global' setting; expected.)
poly_cor_all   <- polychoric(dataallord, correct = FALSE)
cor_matrix_all <- poly_cor_all$rho

# Parallel analysis / scree plot to guide the number of factors
png(file.path(out_dir, "FigS1_scree_parallel_allstudents.png"),
    width = 8, height = 6, units = "in", res = 300)
fa.parallel(cor_matrix_all, cor = "poly", fa = "fa", n.obs = nrow(dataallord))
dev.off()

# 8-factor EFA: principal axis factoring, oblique (oblimin) rotation
fa_results_all <- fa(r = cor_matrix_all,
                     nfactors = 8,
                     n.obs = nrow(dataallord),
                     cor = "poly",
                     rotate = "oblimin",
                     fm = "pa")

print(fa_results_all$loadings, cutoff = 0.3)

loadings_df_all <- as.data.frame(unclass(fa_results_all$loadings))
loadings_df_all$Item <- item_names
loadings_df_all <- relocate(loadings_df_all, Item)
write.csv(loadings_df_all,
          file.path(out_dir, "Supp_EFA_loadings_allstudents.csv"),
          row.names = FALSE)

cat("\n8-factor EFA (all students): RMSEA =", fa_results_all$RMSEA[1],
    " TLI =", fa_results_all$TLI,
    " cumulative variance =", get_cum_var(fa_results_all), "\n")


# ---- 1.2 Citizens / permanent residents: choose number of factors ------------

datacit <- read_data(path_cit)
stopifnot(identical(setdiff(names(datacit), id_cols_cit), item_names))
datacitord <- datacit[, item_names]

poly_cor_cit   <- polychoric(datacitord, correct = FALSE)
cor_matrix_cit <- poly_cor_cit$rho

png(file.path(out_dir, "FigS2_scree_parallel_citizens.png"),
    width = 8, height = 6, units = "in", res = 300)
fa.parallel(cor_matrix_cit, cor = "poly", fa = "fa", n.obs = nrow(datacitord))
dev.off()

# Compare 1-, 2- and 3-factor solutions
fa_fit_compare <- list()
for (k in 1:3) {
  fa_fit_compare[[paste0("f", k)]] <- fa(
    r = cor_matrix_cit,
    nfactors = k,
    n.obs = nrow(datacitord),
    cor = "poly",
    rotate = if (k == 1) "none" else "oblimin",
    fm = "pa"
  )
}

fit_comparison <- data.frame(
  nfactors = 1:3,
  RMSEA    = sapply(fa_fit_compare, function(x) x$RMSEA[1]),
  TLI      = sapply(fa_fit_compare, function(x) x$TLI),
  BIC      = sapply(fa_fit_compare, function(x) x$BIC),
  CumVar   = sapply(fa_fit_compare, get_cum_var)
)
cat("\nCitizens-only EFA: fit comparison across 1-3 factor solutions\n")
print(fit_comparison)
write.csv(fit_comparison,
          file.path(out_dir, "Supp_EFA_fit_comparison_citizens.csv"),
          row.names = FALSE)

# DECISION: use the 3-factor solution.
# Fit improves substantially from 1 to 3 factors (RMSEA 0.148 -> 0.111,
# TLI 0.467 -> 0.699, cumulative variance 33% -> 48%), but all solutions fall
# short of conventional "good fit" cut-offs for this heterogeneous 35-item
# battery. Beyond 3 factors the construct fragments, which defeats the purpose
# of a single composite score.
n_factors_citizens <- 3

fa_results_cit <- fa_fit_compare[[paste0("f", n_factors_citizens)]]
print(fa_results_cit$loadings, cutoff = 0.3)


# ---- 1.3 Item weights for the climate score ----------------------------------
# Each item's weight = its loading on its PRIMARY factor (the factor with the
# largest absolute loading). Items whose primary loading is below 0.30 in
# absolute value get weight 0. This gives one composite score without having
# to pick a single "climate" factor or report three separate sub-scores.

load_mat <- unclass(fa_results_cit$loadings)                      # items x factors

# psych::fa() labels factors PA1, PA2, PA3 and the ORDER of the columns can
# differ between package versions, so the factors are identified by content
# rather than by position. Each anchor item loads most strongly on the factor
# it names:
#   Factor 1 = departmental / relational climate
#   Factor 2 = advisor relationship
#   Factor 3 = research / professional engagement
stopifnot(n_factors_citizens == 3)
anchor_items <- c("Would recommend department",
                  "Rate relationship with my advisor",
                  "How often teach semester course or recitation")
factor_order <- vapply(anchor_items,
                       function(it) which.max(abs(load_mat[item_names == it, ])),
                       integer(1))
stopifnot(!anyDuplicated(factor_order))     # the 3 anchors must map to 3 different factors
load_mat <- load_mat[, factor_order]
colnames(load_mat) <- paste0("Factor", 1:3)

primary_factor  <- apply(abs(load_mat), 1, which.max)             # 1, 2 or 3
primary_loading <- load_mat[cbind(seq_len(nrow(load_mat)), primary_factor)]
primary_loading[abs(primary_loading) < 0.30] <- 0
weights <- as.numeric(primary_loading)
weights[is.na(weights)] <- 0

cat("\nWeights used for climate score:\n")
print(data.frame(Item = item_names, PrimaryFactor = primary_factor,
                 Weight = round(weights, 3)))

# Loadings table (Appendix Table A1)
loadings_cit <- as.data.frame(load_mat)
loadings_cit$Item          <- item_names
loadings_cit$PrimaryFactor <- primary_factor
loadings_cit$WeightUsed    <- round(weights, 3)
loadings_cit <- relocate(loadings_cit, Item)
write.csv(loadings_cit,
          file.path(out_dir, "TableA1_citizens_EFA_loadings.csv"),
          row.names = FALSE)


# ------------------------------------------------------------------------------
# 2. Composite climate score and quartiles (citizens / permanent residents)
# ------------------------------------------------------------------------------

# Missing item responses are replaced with the item mean (keeps N = 1,021)
datacitordimp <- datacitord %>%
  mutate(across(everything(), ~ ifelse(is.na(.), mean(., na.rm = TRUE), .)))
cat("\nRemaining NAs after imputation:", sum(is.na(datacitordimp)), "\n")

# Items 1-13 use a 5-point scale; items 14-35 use a 4-point scale.
five_point_questions <- 1:13
four_point_questions <- 14:35
stopifnot(max(datacitordimp[, five_point_questions]) <= 5,
          max(datacitordimp[, four_point_questions]) <= 4)

# Rescale 5-point items onto the 4-point metric: (x - 1) * 3/4 + 1  (1 -> 1, 5 -> 4)
datacitordimp_normalized <- datacitordimp
datacitordimp_normalized[, five_point_questions] <-
  (datacitordimp[, five_point_questions] - 1) * (3 / 4) + 1

# Climate score = weighted mean of the 35 rescaled items (range 1-4)
datacitordimp_normalized$climatescore <- apply(
  datacitordimp_normalized[, item_names], 1,
  function(x) weighted.mean(x, weights, na.rm = TRUE)
)

# Working dataset: climate score + grouping variables
dat <- datacitordimp_normalized %>%
  mutate(
    inst_type = datacit$InstType,   # PWI / HBCU
    gender    = datacit$Gender,
    urm       = datacit$USURM       # URM / NotURM / Prefer not to answer / Other
  )

# Quartiles of the climate score (Q1 = lowest climate ... Q4 = highest)
quartile_cuts <- quantile(dat$climatescore, probs = seq(0, 1, 0.25), na.rm = TRUE)
cat("\nClimate score quartile cut points:\n")
print(quartile_cuts)

dat <- dat %>%
  mutate(
    climate_quartile = cut(climatescore,
                           breaks = quartile_cuts,
                           labels = c("Q1", "Q2", "Q3", "Q4"),
                           include.lowest = TRUE),
    climate_quartile = factor(climate_quartile,
                              levels = c("Q1", "Q2", "Q3", "Q4"), ordered = TRUE)
  )

cat("\nClimate score summary:\n")
print(summary(dat$climatescore))
cat("\nQuartile sizes:\n")
print(table(dat$climate_quartile))


# ---- Appendix Table A2: item means by climate quartile -----------------------
# (rescaled, imputed items -- the same values that enter the climate score)
quartile_item_means <- dat %>%
  group_by(climate_quartile) %>%
  summarise(across(all_of(item_names), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

quartile_item_means_t <- quartile_item_means %>%
  pivot_longer(-climate_quartile, names_to = "Item", values_to = "Mean") %>%
  pivot_wider(names_from = climate_quartile, values_from = Mean)

write.csv(quartile_item_means_t,
          file.path(out_dir, "TableA2_item_means_by_quartile.csv"),
          row.names = FALSE)


# ---- Figure 1: 15 largest Q4 - Q1 gaps in item means -------------------------
gap_df <- quartile_item_means_t %>%
  mutate(Difference = Q4 - Q1) %>%
  arrange(desc(abs(Difference))) %>%
  slice_head(n = 15)

# Wrap long item text onto several lines for the y-axis
wrap_item_labels <- function(x) {
  vapply(x, function(label) paste(strwrap(label, width = 30), collapse = "\n"),
         character(1))
}

fig1 <- ggplot(gap_df,
               aes(x = Difference, y = fct_reorder(Item, Difference))) +
  geom_col(width = 0.55, fill = "#42677D") +
  geom_text(aes(label = sprintf("%.2f", Difference)),
            hjust = -0.25, size = 3.4, color = "#252525") +
  scale_y_discrete(labels = wrap_item_labels, expand = expansion(add = 0.7)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.13))) +
  labs(x = "Difference in mean item score (Q4 - Q1)", y = NULL) +
  theme_minimal(base_size = 11, base_family = "sans") +
  theme(
    axis.text.y  = element_text(size = 10, hjust = 1, vjust = 0.5,
                                lineheight = 0.95, color = "#252525",
                                margin = margin(r = 12)),
    axis.text.x  = element_text(size = 10, color = "#252525"),
    axis.title.x = element_text(size = 11, hjust = 0.5, margin = margin(t = 12)),
    axis.ticks   = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_line(color = "#E5E5E5", linewidth = 0.3),
    axis.line.x  = element_line(color = "#707070", linewidth = 0.35),
    plot.margin  = margin(t = 12, r = 20, b = 12, l = 12)
  )

ggsave(file.path(out_dir, "Figure1_climate_gaps.png"), fig1,
       width = 9, height = 8, units = "in", dpi = 600, bg = "white")
print(fig1)


# ------------------------------------------------------------------------------
# 3. Who experiences good climate
# ------------------------------------------------------------------------------

# Grouping variables
#   gender3   : Female / Male / "Prefer not to say/Other"
#               ("I prefer not to say", "Other" and missing are combined)
#   urm3      : NotURM / URM / "Prefer not to say/Other"
#               ("Prefer not to answer", "Other" and missing are combined)
#   univ_type : PWI vs. HBCU
# (The "3" in gender3/urm3 refers to the three-level coding; the ART ANOVA
#  below restricts both to two levels.)
dat <- dat %>%
  mutate(
    gender3 = case_when(
      gender %in% c("Male", "Female")                ~ gender,
      gender %in% c("I prefer not to say", "Other")  ~ "Prefer not to say/Other",
      is.na(gender)                                  ~ "Prefer not to say/Other",
      TRUE                                           ~ NA_character_
    ),
    urm3 = case_when(
      urm %in% c("URM", "NotURM")                    ~ urm,
      urm %in% c("Prefer not to answer", "Other")    ~ "Prefer not to say/Other",
      is.na(urm)                                     ~ "Prefer not to say/Other",
      TRUE                                           ~ NA_character_
    ),
    univ_type = inst_type
  )

gender3_levels  <- c("Female", "Male", "Prefer not to say/Other")
urm3_levels     <- c("NotURM", "URM", "Prefer not to say/Other")
univtype_levels <- c("PWI", "HBCU")

who_exp_data <- dat %>%
  filter(!is.na(gender3), !is.na(urm3), !is.na(univ_type)) %>%
  mutate(
    gender3   = factor(gender3, levels = gender3_levels),
    urm3      = factor(urm3, levels = urm3_levels),
    univ_type = factor(univ_type, levels = univtype_levels)
  )

cat("\nWho-experiences-good-climate analytic sample size:", nrow(who_exp_data), "\n")
cat("Excluded for missing gender/URM/institution type:",
    nrow(dat) - nrow(who_exp_data), "\n")


# ---- Table 1: climate quartile by Gender x URM x Institution type ------------
tab_3way <- table(who_exp_data$gender3, who_exp_data$urm3,
                  who_exp_data$univ_type, who_exp_data$climate_quartile)

tab_3way_df <- as.data.frame(tab_3way)
names(tab_3way_df) <- c("Gender", "URM", "UnivType", "ClimateQuartile", "Freq")
write.csv(tab_3way_df,
          file.path(out_dir, "Table1_quartile_by_gender_race_insttype.csv"),
          row.names = FALSE)


# ---- Table 2: chi-square test with standardized residuals --------------------
# One grouping variable made of the observed Gender x URM x Institution-type
# combinations, cross-tabulated with climate quartile.
#
# The chi-square test is DESCRIPTIVE here: it summarizes the overall
# association and, through standardized residuals, flags groups that are over-
# or under-represented in a quartile. The primary inferential test is the ART
# ANOVA below. When more than 20% of cells have an expected count < 5 (as with
# the sparse HBCU cells), the p-value is obtained by Monte Carlo simulation.
# Cramer's V is reported as an overall effect size.
who_exp_data <- who_exp_data %>%
  mutate(group_combo = interaction(gender3, urm3, univ_type,
                                   drop = TRUE, sep = "_"))

tab_combo <- table(who_exp_data$group_combo, who_exp_data$climate_quartile)
cat("\nObserved Gender x URM x Institution-type combinations (group sizes):\n")
print(table(who_exp_data$group_combo))

expected_combo <- suppressWarnings(chisq.test(tab_combo)$expected)
pct_low_combo  <- mean(expected_combo < 5) * 100
cat("\nPercent of cells with expected count < 5:", round(pct_low_combo, 1), "%\n")

if (pct_low_combo > 20) {
  chi_combo <- chisq.test(tab_combo, simulate.p.value = TRUE, B = 10000)
  cat("NOTE: >20% of cells have expected count < 5; using Monte Carlo chi-square.\n")
} else {
  chi_combo <- chisq.test(tab_combo)
}
cramers_combo <- CramerV(tab_combo)

# (df is NA when the p-value is simulated)
cat("\nChi-square (Gender x URM x Institution-type combination x Quartile):\n")
cat("chi2 =", round(chi_combo$statistic, 2),
    " df =", ifelse(is.na(chi_combo$parameter), "NA (Monte Carlo)", chi_combo$parameter),
    " p =", format.pval(chi_combo$p.value, digits = 3),
    " Cramer's V =", round(cramers_combo, 3), "\n")
cat("\nStandardized residuals:\n")
print(round(chi_combo$stdres, 2))

write.csv(as.data.frame(tab_combo),
          file.path(out_dir, "Supp_combo_quartile_crosstab.csv"), row.names = FALSE)
write.csv(round(chi_combo$stdres, 2),
          file.path(out_dir, "Table2_standardized_residuals.csv"))


# ---- Table 3: three-way ART ANOVA --------------------------------------------
# URM x Gender x Institution type (all interactions) on the continuous climate
# score. The full 3 x 3 x 2 design has empty cells once "Prefer not to say/
# Other" is crossed with HBCU, so the 3-way interaction cannot be estimated.
# The ANOVA is therefore run on the 2 x 2 x 2 subsample (Male/Female x
# URM/NotURM x PWI/HBCU), where all 8 cells are populated. Table 1 and the
# chi-square test above still describe the full 3 x 3 x 2 sample.
art_data_2x2x2 <- who_exp_data %>%
  filter(gender3 %in% c("Female", "Male"), urm3 %in% c("NotURM", "URM")) %>%
  mutate(
    gender3   = factor(as.character(gender3), levels = c("Female", "Male")),
    urm3      = factor(as.character(urm3), levels = c("NotURM", "URM")),
    univ_type = factor(as.character(univ_type), levels = c("PWI", "HBCU"))
  )

cat("\n3-way ART ANOVA analytic sample (2x2x2 design, n =",
    nrow(art_data_2x2x2), ")\n")

cell_sizes_3way <- art_data_2x2x2 %>%
  count(univ_type, urm3, gender3, name = "n") %>%
  arrange(univ_type, urm3, gender3)
print(cell_sizes_3way)
write.csv(cell_sizes_3way,
          file.path(out_dir, "Supp_ART_cell_sizes.csv"), row.names = FALSE)

art_model_3way <- art(climatescore ~ urm3 * gender3 * univ_type,
                      data = art_data_2x2x2)
art_anova_3way <- anova(art_model_3way)
cat("\n3-way ART ANOVA: URM x Gender x Institution type on climate score\n")
print(art_anova_3way)
write.csv(as.data.frame(art_anova_3way),
          file.path(out_dir, "Table3_ART_anova.csv"))


# ------------------------------------------------------------------------------
# 4. Numbers reused by the manuscript tables (read by manuscript_tables.Rmd)
# ------------------------------------------------------------------------------

# Sample sizes: written once so the Rmd never hardcodes an N
sample_sizes <- data.frame(
  sample = c("citizens_PR", "who_exp_analytic", "art_2x2x2",
             paste0("quartile_", levels(dat$climate_quartile))),
  n      = c(nrow(dat), nrow(who_exp_data), nrow(art_data_2x2x2),
             as.integer(table(dat$climate_quartile)))
)
write.csv(sample_sizes, file.path(out_dir, "Sample_sizes.csv"), row.names = FALSE)

# Fit statistics and quartile cut points quoted in table footnotes / text
fit3 <- fa_fit_compare[[paste0("f", n_factors_citizens)]]
key_statistics <- data.frame(
  statistic = c("efa_all8_RMSEA", "efa_all8_TLI", "efa_all8_cumvar",
                "efa_cit3_RMSEA", "efa_cit3_TLI", "efa_cit3_cumvar",
                "quartile_cut_25", "quartile_cut_50", "quartile_cut_75",
                "chisq_combo", "cramers_v_combo"),
  value = c(fa_results_all$RMSEA[1], fa_results_all$TLI, get_cum_var(fa_results_all),
            fit3$RMSEA[1], fit3$TLI, get_cum_var(fit3),
            unname(quartile_cuts["25%"]), unname(quartile_cuts["50%"]),
            unname(quartile_cuts["75%"]),
            unname(chi_combo$statistic), cramers_combo)
)
write.csv(key_statistics, file.path(out_dir, "Key_statistics.csv"), row.names = FALSE)

# Record package versions (results depend on them, notably psych/GPArotation)
pkgs <- c("psych", "GPArotation", "dplyr", "tidyr", "ggplot2", "forcats",
          "DescTools", "ARTool")
writeLines(c(R.version.string,
             paste(pkgs, vapply(pkgs, function(p) as.character(packageVersion(p)),
                                character(1)))),
           file.path(out_dir, "package_versions.txt"))

cat("\nDone. Output written to '", out_dir, "/'.\n", sep = "")
