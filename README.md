# 🧠 ELA × Brain Controllability Analyses (ABCD Study)

This repository contains the full R-based analysis pipeline for the manuscript.

> This study used data from the ABCD Study® (N = 7,970) and investigated how early life adversity (ELA) shapes brain network controllability and, in turn, adolescent behavior and cognition. The analysis incorporates mediation, gene–environment interaction (PRS), and longitudinal modeling using cross-lagged panel models.

---

## 📁 Script Overview

| Script | Description |
|--------|-------------|
| `1.merge_controllability_ELAs.R` | Merges brain controllability metrics with ELA factor scores and covariates. |
| `2.Factor_analysis_ELAs_7970subjs.R` | Performs EFA/CFA to derive five latent ELA factors from 67 indicators. |
| `3.merge_lmm_analysis.R` | Runs linear mixed-effects models to test ELA–brain–behavior associations. |
| `4.mediation_by_lavaan.R` | Tests whether controllability mediates ELA's effects on behavior/cognition. |
| `5.moderated_mediation.R` | Conducts moderated mediation (PROCESS model 59 logic) stratified by PRS for ADHD/ASD. |
| `6.CLPM_analysis.R` | Runs longitudinal CLPM to assess brain–behavior coupling over time. |

---

## 📦 Dependencies

All scripts were developed using **R version ≥ 4.2.0**. Key packages used:
- `lavaan`
- `lme4`
- `tidyverse`
- `psych`
- `ggplot2`
- `semTools`

To install dependencies, run:
```r
install.packages(c("lavaan", "lme4", "tidyverse", "psych", "ggplot2", "semTools"))
