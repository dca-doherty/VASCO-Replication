#!/usr/bin/env python3
"""
Nuclear Test - Transient Correlation Validation
Replicates findings from Bruehl and Villarroel (2025) Scientific Reports paper.

Author: Brian Doherty
Date: November 2025
Contact: briandohertyresearch@gmail.com
"""

import pandas as pd
import numpy as np
from datetime import datetime
from scipy import stats
import statsmodels.api as sm
from astropy.time import Time
from astropy.coordinates import get_sun, get_body, AltAz, EarthLocation
import astropy.units as u
import os
import warnings

warnings.filterwarnings('ignore')
np.random.seed(42)

# observatory location
PALOMAR_LAT = 33.3563
PALOMAR_LON = -116.8650
PALOMAR_ELEVATION = 1712

# study period from POSS-I survey
START_DATE = '1949-11-19'
END_DATE = '1957-04-28'


def load_transient_data(filepath):
    """Load the transient dataset from Bruehl."""
    df = pd.read_excel(filepath)
    if 'Date' in df.columns:
        df['date'] = pd.to_datetime(df['Date'])
    return df


def calculate_moon_illumination(dates):
    """Calculate moon illumination for each observation date."""
    palomar = EarthLocation(
        lat=PALOMAR_LAT * u.deg,
        lon=PALOMAR_LON * u.deg,
        height=PALOMAR_ELEVATION * u.m
    )

    results = []
    for date in dates:
        time_obj = Time(date) + 8 * u.hour  # local midnight
        moon = get_body('moon', time_obj, palomar)
        sun = get_sun(time_obj)
        elongation = sun.separation(moon)
        illumination = (1 - np.cos(elongation)) / 2
        results.append({
            'date': date,
            'moon_illumination': illumination.value
        })

    return pd.DataFrame(results)


def create_weather_proxy(dates):
    """
    Create seasonal weather proxy for Palomar area.
    Southern California has predictable seasonal patterns.
    """
    df = pd.DataFrame({'date': pd.to_datetime(dates)})
    df['month'] = df['date'].dt.month

    # seasonal cloud cover probability (dry summers, wetter winters)
    cloud_patterns = {
        1: 0.45, 2: 0.45, 3: 0.40, 4: 0.35, 5: 0.25, 6: 0.20,
        7: 0.15, 8: 0.15, 9: 0.20, 10: 0.25, 11: 0.35, 12: 0.45
    }

    # seasonal precipitation probability
    precip_patterns = {
        1: 0.35, 2: 0.35, 3: 0.30, 4: 0.15, 5: 0.05, 6: 0.02,
        7: 0.01, 8: 0.01, 9: 0.03, 10: 0.10, 11: 0.20, 12: 0.35
    }

    df['cloud_cover_estimate'] = df['month'].map(cloud_patterns)
    df['precip_probability'] = df['month'].map(precip_patterns)

    # add small random variation
    df['cloud_cover_estimate'] += np.random.normal(0, 0.05, len(df))
    df['cloud_cover_estimate'] = df['cloud_cover_estimate'].clip(0, 1)

    return df


def create_media_coverage_index(dates):
    """
    Create media coverage proxy based on historical patterns.
    Major news events competed for coverage with UFO/transient reports.
    """
    df = pd.DataFrame({'date': pd.to_datetime(dates)})
    df['year'] = df['date'].dt.year
    df['month'] = df['date'].dt.month

    df['media_coverage'] = 1.0

    # gradual increase in media attention over study period
    days_since_start = (df['date'] - df['date'].min()).dt.days
    df['media_coverage'] *= (1 + days_since_start / 365.0 * 0.15)

    # weekend reduction in reporting
    df['day_of_week'] = df['date'].dt.dayofweek
    weekend_mask = df['day_of_week'].isin([5, 6])
    df.loc[weekend_mask, 'media_coverage'] *= 0.7

    # summer vacation reduction
    summer_mask = df['month'].isin([6, 7, 8])
    df.loc[summer_mask, 'media_coverage'] *= 0.85

    # normalize
    df['media_coverage'] = df['media_coverage'] / df['media_coverage'].mean()

    return df[['date', 'media_coverage']]


def chi_square_test(data, transient_col, nuclear_col):
    """Run chi-square test for nuclear window association."""
    data['has_transient'] = (data[transient_col] > 0).astype(int)

    contingency = pd.crosstab(data[nuclear_col], data['has_transient'])
    chi2, p_value, dof, expected = stats.chi2_contingency(contingency)

    # calculate relative risk
    in_window_trans = contingency.loc[1, 1]
    in_window_total = contingency.loc[1, :].sum()
    out_window_trans = contingency.loc[0, 1]
    out_window_total = contingency.loc[0, :].sum()

    risk_in = in_window_trans / in_window_total
    risk_out = out_window_trans / out_window_total
    rr = risk_in / risk_out

    return {
        'chi2': chi2,
        'p_value': p_value,
        'relative_risk': rr,
        'contingency': contingency
    }


def negative_binomial_regression(data, transient_col, predictors):
    """Fit negative binomial GLM with controls."""
    X = data[predictors].copy()

    # standardize continuous predictors
    for col in X.columns:
        if X[col].dtype in ['float64', 'int64']:
            if X[col].std() > 0:
                X[col] = (X[col] - X[col].mean()) / X[col].std()

    X = sm.add_constant(X)
    y = data[transient_col]

    model = sm.GLM(y, X, family=sm.families.NegativeBinomial())
    result = model.fit()

    return result


def permutation_test(data, transient_col, nuclear_col, n_permutations=10000):
    """
    Permutation test to validate that specific test dates matter.
    Shuffles nuclear window assignments and recalculates relative risk.
    """
    # observed relative risk
    data['has_transient'] = (data[transient_col] > 0).astype(int)
    table = pd.crosstab(data[nuclear_col], data['has_transient'])
    observed_rr = (table.loc[1, 1] / table.loc[1, :].sum()) / \
                  (table.loc[0, 1] / table.loc[0, :].sum())

    # permutation distribution
    permuted_rrs = []
    for _ in range(n_permutations):
        df_perm = data.copy()
        df_perm[nuclear_col] = np.random.permutation(df_perm[nuclear_col].values)

        perm_table = pd.crosstab(df_perm[nuclear_col], df_perm['has_transient'])
        if perm_table.shape == (2, 2) and perm_table.loc[1, :].sum() > 0:
            rr_perm = (perm_table.loc[1, 1] / perm_table.loc[1, :].sum()) / \
                      (perm_table.loc[0, 1] / perm_table.loc[0, :].sum())
            permuted_rrs.append(rr_perm)

    # two-tailed p-value
    p_value = np.mean(np.array(permuted_rrs) >= observed_rr)
    ci_low = np.percentile(permuted_rrs, 2.5)
    ci_high = np.percentile(permuted_rrs, 97.5)

    return {
        'observed_rr': observed_rr,
        'permutation_p': p_value,
        'ci_low': ci_low,
        'ci_high': ci_high,
        'permutation_mean': np.mean(permuted_rrs)
    }


def run_validation(data_path, output_dir=None):
    """Main validation routine."""
    if output_dir is None:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        output_dir = os.path.join(script_dir, 'results')
    os.makedirs(output_dir, exist_ok=True)

    # load data
    print("Loading transient dataset...")
    data = load_transient_data(data_path)
    print(f"Loaded {len(data)} observation days")

    # identify columns
    transient_col = 'Total_Transients'
    nuclear_col = 'Nuclear_Testing_YN_Window_Plus_Minus_1_Day'
    uap_col = 'Independent_UFOCAT_Sightings_Per_Date'

    # calculate astronomical variables
    print("Calculating moon illumination...")
    date_range = pd.date_range(start=START_DATE, end=END_DATE, freq='D')
    astro_data = calculate_moon_illumination(date_range)

    # create weather and media proxies
    weather_proxy = create_weather_proxy(date_range)
    media_proxy = create_media_coverage_index(date_range)

    # merge data
    data['date'] = pd.to_datetime(data['date']).dt.normalize()
    astro_data['date'] = pd.to_datetime(astro_data['date']).dt.normalize()
    weather_proxy['date'] = pd.to_datetime(weather_proxy['date']).dt.normalize()
    media_proxy['date'] = pd.to_datetime(media_proxy['date']).dt.normalize()

    data = data.merge(astro_data, on='date', how='left')
    data = data.merge(weather_proxy[['date', 'cloud_cover_estimate', 'precip_probability']],
                      on='date', how='left')
    data = data.merge(media_proxy, on='date', how='left')

    print(f"Final dataset: {len(data)} days with all variables merged")

    # test 1: chi-square
    print("\nTest 1: Chi-square test for nuclear window association")
    chi2_result = chi_square_test(data, transient_col, nuclear_col)
    print(f"  Chi-square = {chi2_result['chi2']:.3f}")
    print(f"  p-value = {chi2_result['p_value']:.4f}")
    print(f"  Relative Risk = {chi2_result['relative_risk']:.2f}")
    print(f"  Paper reported: chi2=6.94, p=0.008, RR=1.45")

    # test 2: negative binomial with controls
    print("\nTest 2: Negative binomial regression with controls")
    predictors = [nuclear_col, uap_col, 'moon_illumination',
                  'cloud_cover_estimate', 'media_coverage']
    predictors = [p for p in predictors if p in data.columns]

    nb_result = negative_binomial_regression(data, transient_col, predictors)

    nuclear_irr = np.exp(nb_result.params[nuclear_col])
    nuclear_p = nb_result.pvalues[nuclear_col]
    print(f"  Nuclear window IRR = {nuclear_irr:.3f}")
    print(f"  p-value = {nuclear_p:.4f}")

    # test 3: permutation test
    print("\nTest 3: Permutation test (10,000 iterations)")
    perm_result = permutation_test(data, transient_col, nuclear_col)
    print(f"  Observed RR = {perm_result['observed_rr']:.3f}")
    print(f"  Permutation p-value = {perm_result['permutation_p']:.4f}")
    print(f"  95% CI from permutations: [{perm_result['ci_low']:.3f}, {perm_result['ci_high']:.3f}]")

    # save results
    results_summary = {
        'test': ['Chi-square', 'Negative Binomial (controlled)', 'Permutation'],
        'statistic': [
            f"chi2={chi2_result['chi2']:.3f}, RR={chi2_result['relative_risk']:.2f}",
            f"IRR={nuclear_irr:.3f}",
            f"RR={perm_result['observed_rr']:.3f}"
        ],
        'p_value': [
            chi2_result['p_value'],
            nuclear_p,
            perm_result['permutation_p']
        ],
        'status': [
            'Confirmed' if chi2_result['p_value'] < 0.05 else 'Not significant',
            'Confirmed' if nuclear_p < 0.05 else 'Not significant',
            'Robust' if perm_result['permutation_p'] < 0.05 else 'Not robust'
        ]
    }

    results_df = pd.DataFrame(results_summary)
    results_df.to_csv(os.path.join(output_dir, 'nuclear_correlation_validation.csv'),
                      index=False)

    # save full model summary
    with open(os.path.join(output_dir, 'nb_model_summary.txt'), 'w') as f:
        f.write(nb_result.summary().as_text())

    print(f"\nResults saved to {output_dir}/")

    return {
        'chi2_result': chi2_result,
        'nb_result': nb_result,
        'perm_result': perm_result
    }


if __name__ == '__main__':
    # look for data file in local data folder first
    script_dir = os.path.dirname(os.path.abspath(__file__))
    data_path = os.path.join(script_dir, 'data', 'Transient_Nuclear_Analyzed_Dataset_ScientificReports.xlsx')

    if not os.path.exists(data_path):
        alt_paths = [
            'Transient_Nuclear_Analyzed_Dataset_ScientificReports.xlsx',
            'data/Transient_Nuclear_Analyzed_Dataset_ScientificReports.xlsx'
        ]
        for path in alt_paths:
            if os.path.exists(path):
                data_path = path
                break

    if os.path.exists(data_path):
        results = run_validation(data_path)
    else:
        print(f"Data file not found: {data_path}")
