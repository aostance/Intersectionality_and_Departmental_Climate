# Intersectionality and Departmental Climate Among STEM Doctoral Students

Replication code and blinded data for the manuscript *Intersectionality and Departmental
Climate Among STEM Doctoral Students in an AGEP Alliance*.

The study asks **who experiences good departmental climate** among STEM doctoral students,
using a composite climate score built from a 35-item survey and comparing it across gender,
racial identification, and institution type (predominantly white institutions, PWI, vs. a
historically Black college/university, HBCU).

## Repository contents

| File | Description |
|---|---|
| `climate_analysis.R` | All quantitative analyses: factor analysis, composite climate score, quartiles, chi-square test, three-way ART ANOVA, Figure 1. Writes results to `output/`. |
| `manuscript_tables.Rmd` | Formats the manuscript tables (Word) from the files in `output/`. |
| `data/students_all_blinded.csv` | All survey respondents (N = 2,486): blinded identifiers + 35 survey items. |
| `data/students_citizens_blinded.csv` | U.S. citizens / permanent residents (N = 1,021): blinded identifiers + gender + racial-identification variable + 35 survey items. |
| `output/` | Created when you run the analysis (tables as CSV, figures as PNG). |

## Data

Both files are CSV with one row per respondent, in random order. Blank cells are missing values.

**Identifier columns**

| Column | Files | Values |
|---|---|---|
| `Institution` | both | Arbitrary codes `Inst1`–`Inst4` (not alphabetical or ordered by size) |
| `InstType` | both | `PWI` or `HBCU` |
| `Department` | both | Arbitrary codes `Dept01`–`Dept34`, one per department. Included for completeness; not used in the analyses |
| `Gender` | citizens only | `Female`, `Male`, `I prefer not to say`, `Other`, blank |
| `USURM` | citizens only | `URM`, `NotURM`, `Prefer not to answer`, `Other`, blank. `URM` = Black, Hispanic, and/or Indigenous (BHI); `NotURM` = White and/or Asian (WA) |

**Survey items (35 columns, same order in both files).** The column header is the item text.
Items 1–13 are on a 5-point scale and items 14–35 on a 4-point scale; 1 is always the
lowest/worst response. **Keep the column order unchanged:** the script identifies the
5-point and 4-point items by position.

**Blinding.** Institutions and departments are replaced by arbitrary codes. The key linking
codes to real names is held by the authors and is not part of this repository. Open-ended
(free-text) responses, original department codes, and demographic or program variables that
the analyses do not use have been removed.

## Requirements

R 4.1 or later (the Rmd uses the native pipe `|>`) and these packages:

```r
install.packages(c("psych", "GPArotation", "dplyr", "tidyr", "ggplot2", "forcats",
                   "DescTools", "ARTool", "flextable", "officer", "rmarkdown"))
```

`GPArotation` is **required**: `psych::fa()` needs it for the oblimin rotation. Without it,
`psych` falls back to an unrotated solution and every downstream result changes.

## How to reproduce

1. Open the repository folder as the working directory (or open it as an RStudio project).
2. Run the analysis (a few minutes; it includes parallel analysis and a Monte Carlo test):
   ```r
   source("climate_analysis.R")
   ```
3. Build the formatted tables (writes `manuscript_tables.docx`):
   ```r
   rmarkdown::render("manuscript_tables.Rmd")
   ```

The script sets a random seed (`set.seed(123)`). Parallel analysis and the Monte Carlo
chi-square use random numbers, so keep the order of the analyses unchanged to reproduce the
Monte Carlo p-value exactly. All other results are deterministic.

## Where each manuscript element comes from

| Manuscript element | File in `output/` |
|---|---|
| Figure 1 (largest Q4–Q1 item gaps) | `Figure1_climate_gaps.png` |
| Table 1 (climate quartile by gender × racial identification × institution type) | `Table1_quartile_by_gender_race_insttype.csv` |
| Table 2 (standardized chi-square residuals) | `Table2_standardized_residuals.csv` |
| Table 3 (three-way ART ANOVA) | `Table3_ART_anova.csv` |
| Appendix Table A1 (factor loadings and item weights) | `TableA1_citizens_EFA_loadings.csv` |
| Appendix Table A2 (item means by climate quartile) | `TableA2_item_means_by_quartile.csv` |
| Fit statistics and cut points quoted in the text (EFA RMSEA/TLI/variance, quartile cut points, chi-square, Cramér's V) | `Key_statistics.csv`, `Supp_EFA_fit_comparison_citizens.csv` |
| Sample sizes | `Sample_sizes.csv` |
| Supplementary: all-student EFA loadings, parallel-analysis plots, ART cell sizes, full crosstab for the chi-square | `Supp_*.csv`, `FigS1_*.png`, `FigS2_*.png` |
| Package versions used | `package_versions.txt` |

## Reproducibility notes

- The climate score weights come from an exploratory factor analysis with oblique rotation.
  Factor order and loadings can differ between versions of `psych` and `GPArotation`, so
  `climate_analysis.R` identifies the three factors by item content (not by position) and
  stops with an error if they cannot be told apart. Compare your `output/package_versions.txt`
  with the versions below if results differ.
- Package versions used for the manuscript results: R version 4.4.0 (2024-04-24 ucrt)
psych 2.6.9
GPArotation 2026.8.2
dplyr 1.2.1
tidyr 1.3.2
ggplot2 4.0.3
forcats 1.0.1
DescTools 0.99.60
ARTool 0.11.2


## License and citation


