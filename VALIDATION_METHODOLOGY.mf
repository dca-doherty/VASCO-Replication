# Methodology: Independent Validation of Bruehl and Villarroel Findings

Brian Doherty, Independent Researcher
January 2026
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

1. Bruehl and Villarroel (2025) "Transient detections in historical photographic plates correlate with nuclear weapons testing" - Scientific Reports
2. Villarroel et al. (2025) "Earth shadow deficit analysis of POSS-I transients" - PASP

The validation used the original transient dataset provided by Dr. Stephen Bruehl and applied independent statistical analysis to confirm the reported correlations.

---

## 1. Nuclear Test Window Correlation

### 1.1 Data Sources

Transient Dataset:
- File: Transient_Nuclear_Analyzed_Dataset_ScientificReports.xlsx
- Source: Provided by Dr. Stephen Bruehl
- Period: November 19, 1949 to April 28, 1957 (2,718 days)
- Total transients: 107,862

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
| Center-plate subset | Filtered to transients within 2 degrees of plate center | 2,718 days (22,309 transients) |
| Analyzed | Final dataset for regression and permutation testing | 2,718 days |

No records were excluded. The center-plate filtering reduces transient counts per day but retains all observation days. Environmental variables (precipitation, moon phase) were merged by date. Missing covariate values were rare (less than 1 percent) and handled via listwise deletion in regression models.

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
- Precipitation: NOAA GHCND API, San Diego area stations (1949-1957)
- Moon illumination: Calculated via Astropy ephemeris (accuracy ~1 arcmin)
- Moon altitude: Calculated via Astropy ephemeris

**Seasonal Proxies:**
- Cloud cover: Southern California seasonal patterns (historical daily cloud data unavailable for 1950s)
- Media coverage: Temporal trend proxy for reporting patterns

Note: The nuclear window effect *strengthens* after controlling for these variables (IRR increases from 1.45 crude to 3.53 controlled), indicating confounders were partially masking the true effect rather than creating a spurious correlation.

### 1.5 Negative Binomial Regression

To address overdispersion in count data and control for confounding variables, a negative binomial generalized linear model was fitted.

Model Results (Incidence Rate Ratios):

| Variable | IRR | 95% CI | p-value |
|----------|-----|--------|---------|
| Nuclear Window | 3.527 | 2.799 - 4.446 | < 0.0001 |
| UAP Reports | 1.253 | 1.204 - 1.304 | < 0.0001 |
| Precipitation | 0.370 | 0.329 - 0.418 | < 0.0001 |
| Moon Phase | 0.202 | 0.174 - 0.234 | < 0.0001 |

Key Finding: The nuclear test effect strengthens when controlling for environmental confounders (IRR = 3.53 vs crude RR = 1.45), indicating confounders were partially masking the true effect size.

### 1.6 Center-of-Plate Validation

To rule out edge artifacts from plate scanning, analysis was restricted to transients within 2 degrees of plate center.

Results with Center-Only Data:
- Total transients retained: 22,309 (44% of original)
- Transient-positive days retained: 306 of 310 (98.7%)
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
- Permutation p-value: 0.004
- 95% CI from permutations: 0.708 - 1.319

Interpretation: Only 0.4% of random shuffles produced an effect as large as observed. The specific nuclear test dates matter - this is not explained by temporal autocorrelation.

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
| VASCO (Palomar) | 22,309 | 50 | 0.22% |
| Geometric expectation | - | ~1.4% | - |

Shadow Deficit: Observed rate (0.22%) is significantly below geometric expectation (~1.4%), consistent with original PASP findings.

### 2.3 Pre-Satellite Verification

All shadow transients predate artificial satellites:

| Date | Years Before Sputnik |
|------|---------------------|
| 1949-02-02 | 8.7 years |
| 1951-02-03 | 6.7 years |
| 1953-01-09 | 4.7 years |
| 1956-04-08 | 1.5 years |

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
| Nuclear window (center-only) | - | < 0.0001 | Confirmed |
| Negative binomial (controlled) | IRR = 3.53 | < 0.0001 | Confirmed |
| Permutation test | RR = 1.45 vs mean 1.01 | 0.004 | Confirmed |
| Earth shadow deficit | 0.22% vs ~1.4% expected | Significant | Confirmed |

---

## 5. Technical Appendix

### 5.1 Negative Binomial Model Specification

**Parameterization:** NB2 (quadratic variance function), as implemented in statsmodels.

**Link function:** Log link (canonical for count data).

**Model formula:**
```
log(E[Transients]) = β0 + β1*NuclearWindow + β2*UAP + β3*MoonIllumination 
                     + β4*CloudCover + β5*MediaCoverage
```

**Variable coding:**
- NuclearWindow: Binary (1 = within ±1 day of test, 0 = otherwise)
- UAP: Count (independent UFOCAT sightings per date)
- MoonIllumination: Continuous (0-1 scale from Astropy ephemeris)
- CloudCover: Continuous (0-1 seasonal estimate)
- MediaCoverage: Continuous (normalized temporal trend)

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

**Window definition:** The ±1 day nuclear test window was pre-specified by Bruehl and Villarroel (2025) in the original Scientific Reports paper. This validation adopted the same window definition to test replicability. No window shopping or post-hoc optimization was performed.

**Subgroup analyses:** Center-of-plate filtering (within 2°) was pre-specified by Dr. Villarroel as a methodological refinement to address edge artifact concerns. This was not a data-driven selection.

**Reported p-values:** All p-values reported without correction because each test addresses a distinct hypothesis:
1. Chi-square: Does the association exist?
2. Negative binomial: Does it persist after controlling for confounders?
3. Permutation: Is it driven by specific test dates rather than temporal autocorrelation?

These are complementary validations, not multiple comparisons on the same hypothesis.

### 5.5 Sensitivity Analyses

**Window size:** Original paper tested ±1 day window. Effect persists with ±2 day and ±4 day windows, though diluted as expected (more non-signal days included).

**Temporal subset:** Effect present in both early period (1949-1953) and late period (1953-1957), ruling out single-era artifacts.

**Covariate exclusion:** Nuclear window effect remains significant (p < 0.01) when each covariate is individually removed from the model, confirming robustness.

**Alternative models tested:**
| Model | Nuclear Window Effect | Significant? |
|-------|----------------------|--------------|
| Poisson GLM | IRR = 3.41 | Yes (p < 0.0001) |
| Negative Binomial | IRR = 3.53 | Yes (p < 0.0001) |
| Zero-Inflated Poisson | IRR = 3.38 | Yes (p < 0.0001) |
| Hurdle Model | IRR = 3.47 | Yes (p < 0.0001) |

All model specifications yield consistent estimates, confirming the effect is not an artifact of distributional assumptions.

### 5.6 Earth Shadow Classification Algorithm

**Coordinate inputs:**
- Transient RA/Dec: J2000 coordinates from VASCO catalog
- Observation time: UTC timestamp from plate metadata

**Sun position:** Calculated using Meeus algorithm (Astronomical Algorithms, 2nd ed.) with accuracy ~1 arcminute.

**Anti-sun position:** Shadow center = Sun RA + 180°, Dec = -Sun Dec.

**Shadow geometry:**
- Earth radius: 6,371 km
- GEO orbital radius: 42,164 km (from Earth center)
- Umbra length: 1,380,000 km
- Umbra radius at GEO: 6,234 km
- Shadow angular radius: arcsin(6,234 / 42,164) = 8.50°

**Classification:** Transient classified as "in shadow" if angular separation from anti-sun < 8.50°.

**Validation:** Algorithm tested against JPL Horizons ephemeris for 100 random dates; all shadow center positions agreed within 0.02°.

### 5.7 Limitations

1. **Single observatory:** Palomar-only data cannot rule out site-specific artifacts. Cross-validation with European observatories (Hamburg, in progress) addresses this.

2. **Proxy variables:** Cloud cover and media coverage are seasonal estimates, not daily measurements. However, the nuclear effect *strengthens* when these controls are added, indicating they are not inflating the association.

3. **Transient definition:** This validation uses the transient catalog as provided by the original authors. Any systematic errors in transient detection would propagate to this analysis.

4. **Pre-satellite era:** Shadow analysis assumes no artificial satellites existed during the observation period (1949-1957). This is factually correct (Sputnik launched October 1957), but limits applicability to historical data.

5. **Causal inference:** Statistical association does not establish causation. The finding that transients correlate with nuclear test timing is robust, but the mechanism remains unexplained.

### 5.8 Independence Statement

I conducted this validation on my own time, using my own computing resources. I have no financial interest in the outcome and received no funding or compensation related to this work. I am not affiliated with any university astronomy department, government agency, or organization involved in UAP research.

My only connection to this project is that Dr. Stephen Bruehl shared the transient dataset with me after I reached out expressing interest in independently verifying the published findings. I had no involvement in the original data collection, transient detection pipeline, or the writing of the Scientific Reports or PASP papers. I did not have access to any unpublished data beyond what was shared for validation purposes.

The analysis scripts and output files are available upon request. The underlying transient data belongs to the original research team and should be requested directly from Dr. Bruehl or Dr. Villarroel.

I approached this validation hoping to either confirm or refute the findings. The data confirmed them. I have tried to document the methodology clearly enough that someone else could repeat this work and arrive at the same conclusions, or find errors in my approach if they exist.

---

## 6. File Locations

This validation package is self-contained in this folder.

Validation Scripts:
- nuclear_transient_correlation.py
- earth_shadow_validation.py

Data Files (in data/ subfolder):
- Transient_Nuclear_Analyzed_Dataset_ScientificReports.xlsx
- SUPERVIKTIG_HELAVASCO.csv
- SUPERVIKTIG_HELAVASCO_validated_v4.csv

Output Files (in results/ subfolder after running scripts):
- nuclear_correlation_validation.csv
- nb_model_summary.txt
- shadow_classification.csv
- umbra_transients.csv

---

## 7. Contact

Brian Doherty
Email: briandohertyresearch@gmail.com
Date: January 19, 2026
