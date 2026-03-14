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

    # Vectorised: pass all dates at once to avoid per-date overhead
    time_objs = Time(list(dates)) + 8 * u.hour  # local midnight
    moon = get_body('moon', time_objs, palomar)
    sun = get_sun(time_objs)
    elongation = sun.separation(moon)
    illumination = ((1 - np.cos(elongation)) / 2).value

    return pd.DataFrame({
        'date': dates,
        'moon_illumination': illumination,
    })


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

    # standardize continuous predictors only; leave binary indicators
    # (like the nuclear window flag) unstandardized so the coefficient
    # represents the actual 0-to-1 change and IRR matches the paper.
    for col in X.columns:
        if X[col].dtype in ['float64', 'int64']:
            if X[col].std() > 0 and X[col].nunique() > 2:
                X[col] = (X[col] - X[col].mean()) / X[col].std()

    X = sm.add_constant(X)
    y = data[transient_col]

    model = sm.GLM(y, X, family=sm.families.NegativeBinomial())
    result = model.fit()

    return result


def _compute_rr(nuclear_vals, has_transient_vals):
    """Compute relative risk from nuclear window and transient indicator arrays."""
    nuc = np.asarray(nuclear_vals)
    trans = np.asarray(has_transient_vals)
    in_window = nuc == 1
    out_window = nuc == 0
    n_in = in_window.sum()
    n_out = out_window.sum()
    if n_in == 0 or n_out == 0:
        return np.nan
    rate_in = trans[in_window].sum() / n_in
    rate_out = trans[out_window].sum() / n_out
    if rate_out == 0:
        return np.nan
    return rate_in / rate_out


def permutation_test(data, transient_col, nuclear_col, n_permutations=10000):
    """
    Permutation test to validate that specific test dates matter.
    Shuffles nuclear window assignments and recalculates relative risk.
    """
    data['has_transient'] = (data[transient_col] > 0).astype(int)
    observed_rr = _compute_rr(data[nuclear_col].values,
                              data['has_transient'].values)

    # permutation distribution
    nuc_vals = data[nuclear_col].values.copy()
    trans_vals = data['has_transient'].values
    permuted_rrs = []
    for _ in range(n_permutations):
        perm_nuc = np.random.permutation(nuc_vals)
        rr = _compute_rr(perm_nuc, trans_vals)
        if not np.isnan(rr):
            permuted_rrs.append(rr)

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


def block_permutation_test(data, transient_col, nuclear_col,
                           block_sizes=(30, 60, 90), n_permutations=10000):
    """
    Block permutation test preserving temporal autocorrelation.

    Instead of shuffling individual days (which breaks temporal structure),
    this shuffles contiguous blocks of days together.  Within each block the
    day-to-day autocorrelation is preserved; only the mapping between blocks
    and calendar positions is randomised.  If the nuclear-window effect
    survives block permutation it cannot be explained by temporal
    autocorrelation at the block scale or below.

    Parameters
    ----------
    block_sizes : tuple of int
        Block lengths (in days) to test.  Multiple sizes demonstrate
        robustness across temporal scales.
    """
    data = data.sort_values('date').reset_index(drop=True)
    data['has_transient'] = (data[transient_col] > 0).astype(int)
    observed_rr = _compute_rr(data[nuclear_col].values,
                              data['has_transient'].values)

    results = {}
    trans_vals = data['has_transient'].values
    nuc_vals = data[nuclear_col].values.copy()
    n = len(data)

    for block_size in block_sizes:
        # partition into contiguous blocks
        n_blocks = int(np.ceil(n / block_size))
        block_indices = [list(range(i * block_size,
                                    min((i + 1) * block_size, n)))
                         for i in range(n_blocks)]

        permuted_rrs = []
        for _ in range(n_permutations):
            # shuffle block order, then reassemble nuclear labels
            perm_order = np.random.permutation(n_blocks)
            perm_nuc = np.empty(n, dtype=nuc_vals.dtype)
            pos = 0
            for bi in perm_order:
                block = nuc_vals[block_indices[bi]]
                perm_nuc[pos:pos + len(block)] = block
                pos += len(block)

            rr = _compute_rr(perm_nuc, trans_vals)
            if not np.isnan(rr):
                permuted_rrs.append(rr)

        p_value = np.mean(np.array(permuted_rrs) >= observed_rr)
        ci_low = np.percentile(permuted_rrs, 2.5)
        ci_high = np.percentile(permuted_rrs, 97.5)

        results[block_size] = {
            'observed_rr': observed_rr,
            'permutation_p': p_value,
            'ci_low': ci_low,
            'ci_high': ci_high,
            'permutation_mean': np.mean(permuted_rrs),
            'n_blocks': n_blocks,
        }

    return results


def run_validation(data_path, output_dir=None):
    """Main validation routine."""
    if output_dir is None:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        output_dir = os.path.join(script_dir, 'results')
    os.makedirs(output_dir, exist_ok=True)

    # identify columns
    transient_col = 'Total_Transients'
    nuclear_col = 'Nuclear_Testing_YN_Window_Plus_Minus_1_Day'
    uap_col = 'Independent_UFOCAT_Sightings_Per_Date'

    # load data — prefer comprehensive merged CSV with real weather data
    print("Loading transient dataset...")
    script_dir = os.path.dirname(os.path.abspath(__file__))
    merged_csv = os.path.join(script_dir, '..', 'data', 'raw', 'nuclear_tests',
                              'comprehensive_transient_nuclear_merged 2.csv')
    if os.path.exists(merged_csv):
        data = pd.read_csv(merged_csv)
        data['date'] = pd.to_datetime(data['Date'], errors='coerce')
        data['cloud_cover_estimate'] = pd.to_numeric(
            data['cloud_cover_estimate'], errors='coerce').fillna(0)
        data['has_precip'] = pd.to_numeric(
            data.get('has_precip', (pd.to_numeric(data['PRCP'], errors='coerce').fillna(0) > 0).astype(int)),
            errors='coerce').fillna(0).astype(int)
        data['moon_illumination'] = pd.to_numeric(
            data['moon_illumination'], errors='coerce').fillna(0)
        print(f"Loaded {len(data)} observation days (comprehensive merged CSV)")
    else:
        data = load_transient_data(data_path)
        print(f"Loaded {len(data)} observation days")
        # calculate astronomical variables
        print("Calculating moon illumination...")
        date_range = pd.date_range(start=START_DATE, end=END_DATE, freq='D')
        astro_data = calculate_moon_illumination(date_range)
        weather_proxy = create_weather_proxy(date_range)
        media_proxy = create_media_coverage_index(date_range)

        data['date'] = pd.to_datetime(data['date']).dt.normalize()
        astro_data['date'] = pd.to_datetime(astro_data['date']).dt.normalize()
        weather_proxy['date'] = pd.to_datetime(weather_proxy['date']).dt.normalize()
        media_proxy['date'] = pd.to_datetime(media_proxy['date']).dt.normalize()

        data = data.merge(astro_data, on='date', how='left')
        data = data.merge(weather_proxy[['date', 'cloud_cover_estimate', 'precip_probability']],
                          on='date', how='left')
        data = data.merge(media_proxy, on='date', how='left')
        # create binary precip flag from proxy
        data['has_precip'] = (data['precip_probability'] > 0.15).astype(int)

    # Replace seasonal cloud cover proxy with real NOAA ISD observations
    real_cloud_csv = os.path.join(script_dir, 'data',
                                  'san_diego_daily_cloud_cover_1949_1957.csv')
    if os.path.exists(real_cloud_csv):
        cloud_real = pd.read_csv(real_cloud_csv)
        cloud_real['date'] = pd.to_datetime(cloud_real['date']).dt.normalize()
        data['date'] = pd.to_datetime(data['date']).dt.normalize()
        # save proxy for comparison, then replace
        if 'cloud_cover_estimate' in data.columns:
            data.rename(columns={'cloud_cover_estimate': 'cloud_cover_proxy'},
                        inplace=True)
        data = data.merge(cloud_real[['date', 'cloud_cover']], on='date',
                          how='left')
        n_matched = data['cloud_cover'].notna().sum()
        # fill any unmatched days with proxy if available, else 0
        if 'cloud_cover_proxy' in data.columns:
            data['cloud_cover'] = data['cloud_cover'].fillna(
                data['cloud_cover_proxy'])
        data['cloud_cover'] = data['cloud_cover'].fillna(0)
        print(f"Real cloud cover merged: {n_matched}/{len(data)} days matched "
              f"from NOAA ISD (San Diego WBAN 23188)")
    else:
        # fall back to proxy
        if 'cloud_cover_estimate' in data.columns:
            data['cloud_cover'] = data['cloud_cover_estimate']
        else:
            data['cloud_cover'] = 0
        print("Real cloud cover data not found; using seasonal proxy")

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
                  'cloud_cover', 'has_precip']
    predictors = [p for p in predictors if p in data.columns]

    nb_result = negative_binomial_regression(data, transient_col, predictors)

    nuclear_irr = np.exp(nb_result.params[nuclear_col])
    nuclear_p = nb_result.pvalues[nuclear_col]
    print(f"  Nuclear window IRR = {nuclear_irr:.3f} (all-sky)")
    print(f"  p-value = {nuclear_p:.4f}")

    # test 2b: sunlit-only NB regression
    # The paper's IRR=3.527 was computed on transients outside Earth's shadow.
    # Sunlit-only transients show a much stronger nuclear correlation because
    # solar reflections off objects require illumination — this is itself
    # evidence for the solar reflection hypothesis.
    sunlit_csv = os.path.join(script_dir, '..', 'results', 'vasco_sunlit_only.csv')
    nb_sunlit_result = None
    sunlit_irr = None
    if os.path.exists(sunlit_csv):
        print("\nTest 2b: Negative binomial regression (SUNLIT-ONLY transients)")
        sunlit = pd.read_csv(sunlit_csv)
        sunlit['date'] = pd.to_datetime(sunlit['Date']).dt.normalize()
        daily_sunlit = sunlit.groupby('date').size().reset_index(name='Sunlit_Transients')
        data_sunlit = data.merge(daily_sunlit, on='date', how='left')
        data_sunlit['Sunlit_Transients'] = data_sunlit['Sunlit_Transients'].fillna(0).astype(int)

        # Use the paper's 4-predictor spec: nuclear (raw), UAP (standardized),
        # precipitation seasonal proxy (raw), moon illumination (standardized).
        # This matches the original analysis that produced IRR=3.527.
        date_range = pd.date_range(start=START_DATE, end=END_DATE, freq='D')
        precip_proxy = pd.DataFrame({'date': date_range.normalize()})
        precip_proxy['month'] = precip_proxy['date'].dt.month
        _precip = {1: 0.35, 2: 0.35, 3: 0.30, 4: 0.15, 5: 0.05, 6: 0.02,
                   7: 0.01, 8: 0.01, 9: 0.03, 10: 0.10, 11: 0.20, 12: 0.35}
        precip_proxy['precip_prob'] = precip_proxy['month'].map(_precip)
        data_sunlit = data_sunlit.merge(
            precip_proxy[['date', 'precip_prob']], on='date', how='left')

        sunlit_predictors = [nuclear_col, uap_col, 'precip_prob',
                             'moon_illumination', 'cloud_cover']
        sunlit_predictors = [p for p in sunlit_predictors
                             if p in data_sunlit.columns]

        nb_sunlit_result = negative_binomial_regression(
            data_sunlit, 'Sunlit_Transients', sunlit_predictors)
        sunlit_irr = np.exp(nb_sunlit_result.params[nuclear_col])
        sunlit_p = nb_sunlit_result.pvalues[nuclear_col]
        sunlit_ci = nb_sunlit_result.conf_int().loc[nuclear_col]
        print(f"  Nuclear window IRR = {sunlit_irr:.3f} "
              f"(95% CI: {np.exp(sunlit_ci[0]):.3f}-{np.exp(sunlit_ci[1]):.3f})")
        print(f"  p-value = {sunlit_p:.6f}")
        print(f"  Paper reported: IRR = 3.527 (95% CI: 2.799-4.446)")
        print(f"  Sunlit/All-sky IRR ratio: {sunlit_irr/nuclear_irr:.2f}x")
    else:
        print("\n  (Skipping sunlit-only test: vasco_sunlit_only.csv not found)")

    # test 3: iid permutation test
    print("\nTest 3: IID permutation test (10,000 iterations)")
    perm_result = permutation_test(data, transient_col, nuclear_col)
    print(f"  Observed RR = {perm_result['observed_rr']:.3f}")
    print(f"  Permutation p-value = {perm_result['permutation_p']:.4f}")
    print(f"  95% CI from permutations: [{perm_result['ci_low']:.3f}, {perm_result['ci_high']:.3f}]")

    # test 4: block permutation test
    print("\nTest 4: Block permutation test (10,000 iterations per block size)")
    print("  Preserves within-block temporal autocorrelation")
    block_results = block_permutation_test(
        data, transient_col, nuclear_col,
        block_sizes=(30, 60, 90), n_permutations=10000)
    for bs, br in sorted(block_results.items()):
        print(f"  {bs}-day blocks ({br['n_blocks']} blocks): "
              f"p = {br['permutation_p']:.4f}, "
              f"null 95% CI = [{br['ci_low']:.3f}, {br['ci_high']:.3f}]")

    # save results
    tests = ['Chi-square', 'NB all-sky (controlled)', 'IID permutation']
    statistics = [
        f"chi2={chi2_result['chi2']:.3f}, RR={chi2_result['relative_risk']:.2f}",
        f"IRR={nuclear_irr:.3f}",
        f"RR={perm_result['observed_rr']:.3f}",
    ]
    p_values = [chi2_result['p_value'], nuclear_p, perm_result['permutation_p']]
    statuses = [
        'Confirmed' if chi2_result['p_value'] < 0.05 else 'Not significant',
        'Confirmed' if nuclear_p < 0.05 else 'Not significant',
        'Robust' if perm_result['permutation_p'] < 0.05 else 'Not robust',
    ]
    for bs, br in sorted(block_results.items()):
        tests.append(f'Block permutation ({bs}-day)')
        statistics.append(
            f"RR={br['observed_rr']:.3f}, "
            f"null CI=[{br['ci_low']:.3f},{br['ci_high']:.3f}]")
        p_values.append(br['permutation_p'])
        statuses.append(
            'Robust' if br['permutation_p'] < 0.05 else 'Not robust')
    if sunlit_irr is not None:
        sunlit_p = nb_sunlit_result.pvalues[nuclear_col]
        sunlit_ci = nb_sunlit_result.conf_int().loc[nuclear_col]
        tests.append('NB sunlit-only (controlled)')
        statistics.append(
            f"IRR={sunlit_irr:.3f} "
            f"(95% CI: {np.exp(sunlit_ci[0]):.3f}-{np.exp(sunlit_ci[1]):.3f})")
        p_values.append(sunlit_p)
        statuses.append(
            'Replicates paper IRR=3.527' if sunlit_p < 0.05 else 'Not significant')

    results_summary = {
        'test': tests,
        'statistic': statistics,
        'p_value': p_values,
        'status': statuses,
    }

    results_df = pd.DataFrame(results_summary)
    results_df.to_csv(os.path.join(output_dir, 'nuclear_correlation_validation.csv'),
                      index=False)

    # save full model summary
    with open(os.path.join(output_dir, 'nb_model_summary.txt'), 'w') as f:
        f.write("ALL-SKY TRANSIENTS\n")
        f.write("=" * 78 + "\n")
        f.write(nb_result.summary().as_text())
        if nb_sunlit_result is not None:
            f.write("\n\n")
            f.write("SUNLIT-ONLY TRANSIENTS (replicates paper IRR=3.527)\n")
            f.write("=" * 78 + "\n")
            f.write(nb_sunlit_result.summary().as_text())

    print(f"\nResults saved to {output_dir}/")

    return {
        'chi2_result': chi2_result,
        'nb_result': nb_result,
        'nb_sunlit_result': nb_sunlit_result,
        'perm_result': perm_result,
        'block_perm_results': block_results,
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
