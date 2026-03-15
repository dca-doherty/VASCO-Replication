# Methodology: Independent Validation of Bruehl and Villarroel Findings

Brian Doherty, Independent Researcher
Contact: briandohertyresearch@gmail.com

---

## About the Analyst

I am a Senior Data Analyst with over a decade of experience in statistical modeling, pattern recognition, and large-scale data validation. My professional work involves analyzing complex datasets to identify anomalies, validate business hypotheses, and build predictive models using regression techniques, time-series analysis, and machine learning methods. I hold a degree in Economics from Washington State University.

While my day job is in the financial services industry rather than astronomy, the core analytical techniques are identical: cleaning messy data, controlling for confounding variables, testing for statistical significance, and remaining skeptical of apparent patterns until they survive rigorous validation. I routinely work with datasets containing millions of records, build negative binomial and logistic regression models, conduct permutation testing, and present findings to stakeholders who demand methodological rigor.

This validation was conducted independently in my personal time. I approached the Bruehl and Villarroel findings as I would any anomaly detection problem: assume the null hypothesis, control for known confounders, test whether the signal survives multiple analytical approaches, and document the methodology transparently. The fact that the subject matter involves astronomical transients rather than financial transactions does not change the underlying statistical principles.

I have no institutional affiliation with any astronomy department or UFO/UAP research organization. My only stake in this analysis is intellectual honesty.

---

## Overview

This document describes the methodology used to independently validate findings from two published papers:

1. Bruehl and Villarroel (2025) "Transients in the Palomar Observatory Sky Survey (POSS-I) may be associated with nuclear testing and reports of unidentified anomalous phenomena" - Scientific Reports, 15, 34125
2. Villarroel et al. (2025) "Aligned, Multiple-transient Events in the First Palomar Sky Survey" - PASP, 137, 104504

The validation used the original transient dataset provided by Dr. Stephen Bruehl and applied independent statistical analysis to confirm the reported correlations.

---

## 1. Nuclear Test Window Correlation

### 1.1 Data Sources

Transient Dataset:
- File: Transient_Nuclear_Analyzed_Dataset_ScientificReports.xlsx
- Source: Provided by Dr. Stephen Bruehl
- Period: November 19, 1949 to April 28, 1957 (2,718 days)
- Total transients: 107,875

Nuclear Test Data:
- Source: DOE/NV-209 Rev 16 (US Department of Energy official record)
- Cross-reference: Johnston's Archive
- Window definition: Plus or minus 1 day from detonation date

### 1.2 Data Flow

The transient dataset was received as a preprocessed Excel file from Dr. Stephen Bruehl. I did not perform raw plate extraction or transient detection - that work was done by the VASCO team using their published pipeline (Solano et al. 2022, MNRAS 515, 1380-1391).

| Stage | Description | Records |
|-------|-------------|---------|
| Received | Daily transient counts from Bruehl | 2,718 days |
| Merged | Added nuclear test flags, environmental covariates | 2,718 days |
| Full catalog | All transients classified by shadow position | 107,875 transients |
| Center-plate subset | Filtered to transients within 2 degrees of plate center | 31,525 transients (29.2%) |
| Analyzed | Final dataset for regression and permutation testing | 2,718 days |

No records were excluded from the daily regression. The center-plate filtering reduces transient counts per day but retains all observation days. Plate centers are computed per-plate via unit-vector averaging of source positions, handling RA wraparound at 0/360 degrees. Environmental variables (precipitation, moon phase) were merged by date. Missing covariate values were rare (less than 1 percent) and handled via listwise deletion in regression models.

### 1.3 Chi-Square Analysis

Method: 2x2 contingency table comparing transient detection rates inside versus outside nuclear test windows.

Results:

| Metric | Validation | Original Paper |
|--------|------------|----------------|
| Chi-square | 6.47 | 6.94 |
| p-value | 0.011 | 0.008 |
| Relative Risk | 1.45 | 1.45 |

Contingency Table:

|                    | No Transient | Has Transient | Total | Rate |
|--------------------|--------------|---------------|-------|------|
| Outside Window     | 2,116        | 255           | 2,371 | 10.8% |
| In Nuclear Window  | 293          | 54            | 347   | 15.6% |

Status: CONFIRMED. Relative risk matches exactly at 1.45.

### 1.4 Environmental Control Variables

To rule out confounding by observing conditions, several environmental variables were incorporated into the analysis.

**Real Historical Data:**
- Precipitation: NOAA GHCND API, San Diego area stations (1949-1957), coded as binary (has_precip)
- Cloud cover: Daily mean from NOAA Integrated Surface Database (ISD), San Diego station WBAN 23188, parsed from hourly sky condition reports (GF1 field, oktas converted to 0-1 scale). 3,287 daily records with zero missing days.
- Moon illumination: Calculated via vectorized Astropy ephemeris (accuracy ~1 arcmin)

**Standardization:** Continuous predictors are standardized (zero mean, unit variance). Binary predictors (nuclear window, has_precip) are kept on their natural 0/1 scale so that exp(coefficient) gives the IRR for a 0-to-1 shift.

Note: The nuclear window effect *strengthens* after controlling for these variables (IRR increases from 1.45 crude to 1.80 all-sky controlled), indicating confounders were partially masking the true effect rather than creating a spurious correlation.

### 1.5 Negative Binomial Regression

To address overdispersion in count data and control for confounding variables, a negative binomial generalized linear model was fitted.

**All-sky model results (Incidence Rate Ratios):**

| Variable | IRR | 95% CI | p-value |
|----------|-----|--------|---------|
| Nuclear Window | 1.803 | 1.601 - 2.030 | < 0.0001 |
| UAP Reports | 1.253 | 1.206 - 1.304 | < 0.0001 |
| Precipitation (binary) | 0.369 | 0.327 - 0.416 | < 0.0001 |
| Moon Phase | 0.201 | 0.193 - 0.209 | < 0.0001 |

**Sunlit-only model results (transients outside Earth's shadow):**

| Variable | IRR | 95% CI | p-value |
|----------|-----|--------|---------|
| Nuclear Window | 3.981 | 3.475 - 4.562 | < 0.0001 |
| Original Paper | 3.527 | 2.799 - 4.446 | < 0.0001 |

Key Findings:

The all-sky nuclear window IRR of 1.80 represents an 80% increase in transient counts during nuclear test windows after controlling for environmental confounders.

Restricting to sunlit transients (those outside Earth's geometric shadow) yields IRR = 3.98 (95% CI: 3.475-4.562), which reproduces the original paper's value of 3.53 (95% CI: 2.799-4.446) with overlapping confidence intervals. The original paper's IRR was almost certainly computed on a sunlit-only or shadow-excluded subset.

The near-doubling of the IRR when restricting to sunlit positions (1.80 all-sky versus 3.98 sunlit-only) is itself a physically meaningful result: the nuclear test correlation is concentrated among transients that require solar illumination, independently supporting the solar reflection hypothesis and corroborating the shadow deficit finding.

### 1.6 Center-of-Plate Validation

To rule out edge artifacts from plate scanning, analysis was restricted to transients within 2 degrees of plate center. Plate centers are computed per-plate via unit-vector averaging of source positions.

Results with Center-Only Data:
- Total transients retained: 31,525 (29.2% of original 107,875)
- Nuclear window effect: p < 0.0001

Status: CONFIRMED. Core finding survives edge artifact removal with increased statistical significance.

### 1.7 Permutation Testing

To validate that specific nuclear test dates drive the correlation rather than general temporal autocorrelation:

Procedure:
1. Shuffle nuclear test date assignments randomly
2. Recalculate relative risk for each permutation
3. Compare observed RR to permutation distribution
4. Repeat 10,000 times

Results:
- Observed RR: 1.447
- Permutation mean: 1.009
- Permutation p-value: 0.006
- 95% CI from permutations: 0.708 - 1.319

Interpretation: Only 0.6% of random shuffles produced an effect as large as observed. The specific nuclear test dates matter - this is not explained by temporal autocorrelation.

---

## 2. Earth Shadow Deficit Validation

### 2.1 Methodology

Objects inside Earth's geometric shadow cannot reflect sunlight. If transients were sunlit orbital debris, they should not appear in the shadow region. A deficit of transients in Earth's shadow supports an interpretation inconsistent with debris.

Shadow Parameters:
- Earth radius: 6,371 km
- GEO orbital radius: 42,164 km
- Umbra length: 1,380,000 km
- Umbra radius at GEO: 6,234 km
- Shadow angle (proper umbra): 8.50 degrees
- Shadow angle (simple geometric): 8.69 degrees

Calculation Method:
1. For each transient, compute observation UTC time
2. Calculate Sun position using Meeus algorithm (accuracy ~1 arcmin)
3. Compute anti-sun position (shadow center = Sun RA + 180 deg, Dec = -Sun Dec)
4. Measure angular distance from transient to shadow center
5. Flag transients within shadow radius as "in shadow"

### 2.2 Results

| Dataset | Total | In Shadow | Rate |
|---------|-------|-----------|------|
| VASCO full catalog | 107,875 | 499 | 0.46% |
| VASCO center-of-plate | 31,525 | 142 | 0.45% |
| Geometric expectation | - | - | ~1.4% |

Shadow Deficit: The shadow deficit is present in both the full catalog and the center-of-plate subset, with nearly identical rates (0.46% and 0.45%). The consistency across both filters indicates that edge artifacts do not preferentially populate the shadow region. Both rates fall well below the geometric expectation of approximately 1.4%, consistent with the original PASP findings.

### 2.3 Pre-Satellite Verification

All shadow transients are distributed across nine POSS-I plates and predate artificial satellites:

| Plate | Date | Years Before Sputnik |
|-------|------|---------------------|
| XE472 | 1949-11-19 | 7.9 years |
| XE358 | 1950-12-10 | 6.8 years |
| XE369 | 1951-02-03 | 6.7 years |
| XE356 | 1951-11-29 | 5.9 years |
| XE413 | 1951-11-30 | 5.9 years |
| XE422 | 1953-01-09 | 4.7 years |
| XE355 | 1954-11-22 | 2.9 years |
| XE415 | 1955-12-10 | 1.8 years |
| XE617 | 1956-04-08 | 1.5 years |

First artificial satellite: Sputnik 1, October 4, 1957

Conclusion: Shadow transients cannot be explained by satellite debris since no artificial objects existed during the observation period.

---

## 3. Software Environment

Python Version: 3.12/3.13

Dependencies:
- pandas: Data manipulation
- numpy: Numerical operations
- scipy: Statistical tests
- statsmodels: Regression models
- astropy: Astronomical calculations

Random Seeds: Set to 42 for reproducibility of permutation tests.

---

## 4. Summary of Validated Findings

| Finding | Statistic | p-value | Status |
|---------|-----------|---------|--------|
| Nuclear window chi-square | chi2 = 6.47, RR = 1.45 | 0.011 | Confirmed |
| NB all-sky (controlled) | IRR = 1.80 | < 0.0001 | Confirmed |
| NB sunlit-only (controlled) | IRR = 3.98 (paper: 3.53) | < 0.0001 | Replicates paper |
| Permutation test | RR = 1.45 vs mean 1.01 | 0.006 | Confirmed |
| Earth shadow deficit (full) | 0.46% vs ~1.4% expected | Significant | Confirmed |
| Earth shadow deficit (center) | 0.45% vs ~1.4% expected | Significant | Confirmed |

---

## 5. Technical Appendix

### 5.1 Negative Binomial Model Specification

**Parameterization:** NB2 (quadratic variance function), as implemented in statsmodels.

**Link function:** Log link (canonical for count data).

**All-sky model formula:**

```
log(E[Transients]) = B0 + B1*NuclearWindow + B2*UAP + B3*MoonIllumination
                     + B4*CloudCover + B5*HasPrecip
```

**Sunlit-only model formula:**

```
log(E[SunlitTransients]) = B0 + B1*NuclearWindow + B2*UAP + B3*PrecipProb
                           + B4*MoonIllumination
```

**Variable coding:**
- NuclearWindow: Binary (1 = within +/-1 day of test, 0 = otherwise). Not standardized.
- HasPrecip: Binary (1 = precipitation recorded at NOAA stations, 0 = otherwise). Not standardized.
- UAP: Count (independent UFOCAT sightings per date). Standardized.
- MoonIllumination: Continuous (0-1 scale from Astropy ephemeris). Standardized.
- CloudCover: Continuous (0-1 daily mean from NOAA ISD hourly reports). Standardized.

**Standardization note:** Binary predictors are kept on their natural 0/1 scale so that exp(coefficient) gives the true IRR for being in vs out of that condition. Continuous predictors are standardized (zero mean, unit variance) for numerical stability. Earlier iterations of the script incorrectly standardized all predictors including binary ones, which compressed the nuclear window coefficient.

**No offset term:** Daily counts are already normalized by observation day; no exposure adjustment needed.

**No random effects:** Single observatory (Palomar) eliminates need for station-level clustering.

### 5.2 Model Diagnostics

**Overdispersion test:** Poisson model rejected in favor of negative binomial.
- Poisson deviance/df ratio: 847.3 (severe overdispersion)
- NB2 deviance/df ratio: 1.12 (acceptable)
- Likelihood ratio test (Poisson vs NB): p < 0.0001

**Model comparison (AIC):**

| Model | AIC | Notes |
|-------|-----|-------|
| Poisson | 18,247 | Overdispersed, poor fit |
| Negative Binomial | 12,891 | Preferred |
| Zero-Inflated NB | 12,903 | No improvement; excess zeros not a problem |

**Multicollinearity:** Variance Inflation Factors (VIF) all < 2.0. No problematic collinearity among predictors.

**Influential observations:** Cook's distance computed; no single observation exceeds 0.5. November 10, 1955 (4,528 transients) has highest leverage but removal does not qualitatively change results.

### 5.3 Autocorrelation Assessment

**Temporal autocorrelation:**
- Durbin-Watson statistic on Pearson residuals: 1.87 (acceptable; values near 2.0 indicate no autocorrelation)
- ACF plot shows no significant lags beyond lag-0
- The permutation test explicitly addresses temporal structure by shuffling nuclear test dates while preserving transient time series

**Spatial autocorrelation:** Not applicable. Single observatory (Palomar) means all observations share the same spatial location. Geographic effects tested separately via Hamburg Observatory cross-validation (not included in this validation document).

### 5.4 Multiple Comparisons and Pre-specification

**Window definition:** The +/-1 day nuclear test window was pre-specified by Bruehl and Villarroel (2025) in the original Scientific Reports paper. This validation adopted the same window definition to test replicability. No window shopping or post-hoc optimization was performed.

**Subgroup analyses:** Center-of-plate filtering (within 2 degrees) was pre-specified by Dr. Villarroel as a methodological refinement to address edge artifact concerns. This was not a data-driven selection.

**Reported p-values:** All p-values reported without correction because each test addresses a distinct hypothesis:
1. Chi-square: Does the association exist?
2. Negative binomial: Does it persist after controlling for confounders?
3. Permutation: Is it driven by specific test dates rather than temporal autocorrelation?

These are complementary validations, not multiple comparisons on the same hypothesis.

### 5.5 Sensitivity Analyses

**Window size:** Original paper tested +/-1 day window. Effect persists with +/-2 day and +/-4 day windows, though diluted as expected (more non-signal days included).

**Temporal subset:** Effect present in both early period (1949-1953) and late period (1953-1957), ruling out single-era artifacts.

**Covariate exclusion:** Nuclear window effect remains significant (p < 0.01) when each covariate is individually removed from the model, confirming robustness.

**All-sky vs sunlit-only:** The nuclear window IRR approximately doubles when restricting to sunlit transients (1.80 all-sky versus 3.98 sunlit-only), consistent across all model specifications tested (Poisson, NB1, NB2, zero-inflated Poisson, hurdle model). All specifications yield significant nuclear window effects (p < 0.0001), with all-sky IRR estimates in the range of 1.7-2.0.

### 5.6 Earth Shadow Classification Algorithm

**Coordinate inputs:**
- Transient RA/Dec: J2000 coordinates from VASCO catalog
- Observation time: UTC timestamp from plate metadata

**Sun position:** Calculated using Meeus algorithm (Astronomical Algorithms, 2nd ed.) with accuracy ~1 arcminute.

**Anti-sun position:** Shadow center = Sun RA + 180 degrees, Dec = -Sun Dec.

**Shadow geometry:**
- Earth radius: 6,371 km
- GEO orbital radius: 42,164 km (from Earth center)
- Umbra length: 1,380,000 km
- Umbra radius at GEO: 6,234 km
- Shadow angular radius: arcsin(6,234 / 42,164) = 8.50 degrees

**Classification:** Transient classified as "in shadow" if angular separation from anti-sun < 8.50 degrees.

**Center-of-plate filtering:** Plate centers are computed per-plate via unit-vector averaging of source positions (handles RA wraparound at 0/360 degrees). Transients beyond a configurable radius (default 2.0 degrees) from their plate center are excluded.

**Validation:** Algorithm tested against JPL Horizons ephemeris for 100 random dates; all shadow center positions agreed within 0.02 degrees.

### 5.7 Limitations

1. **Single observatory:** Palomar-only data cannot rule out site-specific artifacts. Cross-validation with European observatories (Hamburg, in progress) addresses this.
2. **Environmental controls:** Cloud cover and precipitation use real NOAA station data from San Diego, approximately 100 km from Palomar. Conditions at the observatory may differ from the coastal station, particularly for marine layer cloud cover. However, the nuclear effect persists with or without these controls, indicating the finding is not driven by environmental confounding.
3. **Transient definition:** This validation uses the transient catalog as provided by the original authors. Any systematic errors in transient detection would propagate to this analysis.
4. **Pre-satellite era:** Shadow analysis assumes no artificial satellites existed during the observation period (1949-1957). This is factually correct (Sputnik launched October 1957), but limits applicability to historical data.
5. **Causal inference:** Statistical association does not establish causation. The finding that transients correlate with nuclear test timing is robust, but the mechanism remains unexplained.

### 5.8 Independence Statement

I conducted this validation on my own time, using my own computing resources. I have no financial interest in the outcome and received no funding or compensation related to this work. I am not affiliated with any university astronomy department, government agency, or organization involved in UAP research.

My only connection to this project is that Dr. Stephen Bruehl shared the transient dataset with me after I reached out expressing interest in independently verifying the published findings. I had no involvement in the original data collection, transient detection pipeline, or the writing of the Scientific Reports or PASP papers. I did not have access to any unpublished data beyond what was shared for validation purposes.

The analysis scripts and output files are publicly available at https://github.com/dca-doherty/VASCO-Replication. The underlying transient data belongs to the original research team and should be requested directly from Dr. Bruehl or Dr. Villarroel.

I approached this validation hoping to either confirm or refute the findings. The data confirmed them. I have tried to document the methodology clearly enough that someone else could repeat this work and arrive at the same conclusions, or find errors in my approach if they exist.

---

## 6. File Locations

This validation package is self-contained in this folder.

Validation Scripts
- nuclear_transient_correlation.py
- earth_shadow_validation.py

Data Files (in data/ subfolder, not included in public repo)
- Transient_Nuclear_Analyzed_Dataset_ScientificReports.xlsx
- SUPERVIKTIG_HELAVASCO.csv
- SUPERVIKTIG_HELAVASCO_validated_v4.csv

Output Files (in results/ subfolder after running)
- nuclear_correlation_validation.csv
- nb_model_summary.txt (all-sky and sunlit-only model output)
- umbra_transients_full.csv (499 rows)
- umbra_transients_center.csv (142 rows, requires --center-plate flag)
