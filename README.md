# Validation Replication Package

Independent validation of core findings from Bruehl and Villarroel (2025), "Transients in the Palomar Observatory Sky Survey (POSS-I) may be associated with nuclear testing and reports of unidentified anomalous phenomena," published in Scientific Reports.

Author: Brian Doherty
Contact: briandohertyresearch@gmail.com
Date: January 2026

--

## Purpose

This package reproduces the primary statistical claims from the Bruehl and Villarroel paper using the original dataset provided by Dr. Stephen Bruehl. The goal is straightforward: take the same data, run independent analyses, and confirm (or fail to confirm) the reported results. Everything here is written from scratch rather than reusing the original authors' code.

Two analyses are included

1. Nuclear test window correlation (chi-square, negative binomial regression, permutation testing)
2. Earth shadow deficit analysis (geometric shadow classification of transient positions)

### Data (Not Included)

The data folder contains the original research datasets and is not included in the public repository per agreement with the data providers. If you need access to the datasets to reproduce the analysis, please contact Dr. Stephen Bruehl at Vanderbilt University Medical Center or Dr. Beatriz Villarroel at Stockholm University.

The scripts expect three files in the data folder:

`Transient_Nuclear_Analyzed_Dataset_ScientificReports.xlsx` Original transient dataset from Dr. Bruehl, containing daily transient counts from POSS-I plates (November 1949 through April 1957, 2,718 days)
`SUPERVIKTIG_HELAVASCO.csv` VASCO transient catalog with source positions and plate metadata
`SUPERVIKTIG_HELAVASCO_validated_v4.csv` Validated subset of the VASCO catalog

## Requirements

Python 3.10 or higher.

Dependencies

pandas
numpy
scipy
statsmodels
astropy
openpyxl (for reading the Excel dataset)

Install with

```bash
pip install pandas numpy scipy statsmodels astropy openpyxl
```

## How to Run

From this directory

```bash
python nuclear_transient_correlation.py
python earth_shadow_validation.py
```

Both scripts print results to the console and save output files to the `results/` subfolder. The random seed is fixed at 42 for reproducibility of permutation tests.

## What the Scripts Do

### nuclear_transient_correlation.py

Loads the original Bruehl dataset and runs three independent tests of the nuclear window correlation:

**Test 1: Chi-square contingency test.** Builds a 2x2 table comparing transient detection rates inside versus outside nuclear test windows (plus or minus 1 day from detonation). Computes chi-square statistic, p-value, and relative risk.

**Test 2: Negative binomial regression.** Fits a generalized linear model with nuclear test window as the primary predictor, controlling for precipitation, cloud cover, moon illumination, and a media coverage proxy. Reports incidence rate ratios with 95% confidence intervals.

**Test 3: Permutation test.** Randomly shuffles nuclear test date assignments 10,000 times, recalculates relative risk each time, and compares the observed value to the null distribution. This confirms whether the specific test dates matter or whether any random set of dates would produce a similar effect.

Output files
`results/nuclear_correlation_validation.csv` -- Summary statistics
`results/nb_model_summary.txt` -- Full negative binomial model output

### earth_shadow_validation.py

Classifies each VASCO transient by whether its sky position falls inside Earth's geometric shadow at geostationary orbit altitude. Objects inside the shadow cannot reflect sunlight, so a deficit of transients in the shadow region has implications for what the transients physically are.

The script calculates the Sun position for each observation time using standard astronomical algorithms, computes the anti-sun point (shadow center), and measures the angular distance from each transient to the shadow center. Transients within the shadow radius at GEO altitude (approximately 8.5 degrees) are flagged.

Output files
`results/shadow_classification.csv` -- Full classification of all transients
`results/umbra_transients.csv` -- Only the transients falling inside the shadow

## Replicated Results

### Nuclear Window Correlation

| Metric | This Validation | Original Paper |
|--------|----------------|----------------|
| Chi-square | 6.47 | 6.94 |
| p-value | 0.011 | 0.008 |
| Relative Risk | 1.45 | 1.45 |

The relative risk matches exactly. The small difference in chi-square likely reflects minor differences in how the contingency table handles edge cases.

### Negative Binomial Regression (with environmental controls)

| Variable | IRR | 95% CI | p-value |
|----------|-----|--------|---------|
| Nuclear Window | 3.527 | 2.799 - 4.446 | < 0.0001 |
| UAP Reports | 1.253 | 1.204 - 1.304 | < 0.0001 |
| Precipitation | 0.370 | 0.329 - 0.418 | < 0.0001 |
| Moon Phase | 0.202 | 0.174 - 0.234 | < 0.0001 |

The nuclear test effect strengthens after controlling for environmental confounders (IRR = 3.53 versus the crude relative risk of 1.45), indicating that weather and lunar conditions were partially masking the true effect size.

### Permutation Test

Observed RR: 1.447
Permutation mean: 1.009
Permutation p-value: 0.004
95% CI from null distribution: 0.708 - 1.319

Only 0.4% of random shuffles produced a relative risk as large as the observed value. The specific nuclear test dates drive the correlation, not temporal autocorrelation.

### Center-of-Plate Validation

Restricting to transients within 2 degrees of plate center (eliminating edge artifacts):
Transients retained: 22,309 (44% of original)
Transient-positive days retained: 306 of 310 (98.7%)
Nuclear window effect: p < 0.0001

The finding survives edge artifact removal with increased statistical significance.

### Earth Shadow Deficit

| Dataset | Total | In Shadow | Rate |
|---------|-------|-----------|------|
| VASCO (Palomar) | 22,309 | 50 | 0.22% |
| Geometric expectation | -- | -- | ~1.4% |

Shadow transient rate (0.22%) is significantly below geometric expectation (~1.4%), consistent with the original findings. All 50 shadow transients predate artificial satellites (earliest: February 1949, latest: April 1956, all before Sputnik in October 1957).

## Methodology

See `VALIDATION_METHODOLOGY.md` in this folder for the complete methodology documentation, including detailed descriptions of the statistical tests, assumptions, shadow geometry calculations, and how results compare to the original paper at each step.

## Data Sources

Nuclear test dates: DOE/NV-209 Rev 16, Johnston's Archive
Transient data: Provided by Dr. Stephen Bruehl (Vanderbilt University Medical Center)
VASCO catalog: VASCO project (Villarroel et al.)
Solar position calculations: Astropy library using standard Meeus algorithms

## Citation

If you use this validation package, please cite:

Bruehl, S. and Villarroel, B. (2025). "Transients in the Palomar Observatory Sky Survey (POSS-I) may be associated with nuclear testing and reports of unidentified anomalous phenomena." Scientific Reports, 15, 34125.

## Contact

Brian Doherty
briandohertyresearch@gmail.com
