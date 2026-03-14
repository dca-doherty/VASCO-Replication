#!/usr/bin/env python3
"""
Plate-by-plate shadow coverage calculation using VASCO 2-degree center-of-plate data.

Uses the SUPERVIKTIG_HELAVASCO_within2deg_CENTER.csv dataset (~22k transients)
which has been pre-filtered to transients within 2 degrees of their plate center.

For each plate, computes:
  1. Plate center (mean RA/Dec of transients on the plate)
  2. Observation time from plate metadata
  3. Anti-sun position at that time
  4. Fraction of the plate's 2-degree circular field overlapping the shadow cone
  5. Sum to get total expected shadow coverage as a fraction of total plate area

Also counts observed transients falling within the shadow cone to compare
expected vs observed shadow transient rates.

Author: Brian Doherty
Date: March 2026
"""

import csv
import math
import os
import sys

# Add parent module
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from earth_shadow_validation import (
    sun_position_meeus, anti_sun_position, angular_separation,
    compute_plate_center, GEO_SHADOW_RADIUS_DEG, DEFAULT_CENTER_PLATE_RADIUS_DEG,
)

PLATE_FIELD_RADIUS_DEG = DEFAULT_CENTER_PLATE_RADIUS_DEG  # 2.0 deg


def circle_overlap_fraction(d, r1, r2):
    """
    Fraction of a circle of radius r1 that overlaps with a circle of radius r2
    when their centers are separated by angular distance d (all in degrees).

    Returns a value between 0 and 1 representing the fraction of the r1-circle
    area that falls inside the r2-circle.
    """
    if d >= r1 + r2:
        return 0.0
    if d + r1 <= r2:
        # plate entirely inside shadow
        return 1.0
    if d + r2 <= r1:
        # shadow entirely inside plate
        return (r2 / r1) ** 2

    # partial overlap: compute area of intersection of two circles
    cos_arg1 = (d * d + r1 * r1 - r2 * r2) / (2 * d * r1)
    cos_arg2 = (d * d + r2 * r2 - r1 * r1) / (2 * d * r2)
    cos_arg1 = max(-1, min(1, cos_arg1))
    cos_arg2 = max(-1, min(1, cos_arg2))

    alpha = math.acos(cos_arg1)
    beta = math.acos(cos_arg2)

    area_intersection = (r1 * r1 * (alpha - math.sin(alpha) * math.cos(alpha))
                         + r2 * r2 * (beta - math.sin(beta) * math.cos(beta)))
    area_plate = math.pi * r1 * r1

    return area_intersection / area_plate


def load_plates(vasco_path):
    """Load VASCO transients and group by plate."""
    plates = {}
    with open(vasco_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            name = row['Name']
            try:
                ra = float(row['RA'])
                dec = float(row['Dec'])
                time_str = row['UTobservation']
            except (ValueError, KeyError):
                continue
            plates.setdefault(name, []).append({
                'ra': ra, 'dec': dec, 'name': name, 'time_str': time_str
            })
    return plates


def parse_obs_time(time_str):
    """Parse observation time string to (year, month, day, hour, min, sec)."""
    from datetime import datetime
    if 'T' in time_str:
        dt = datetime.fromisoformat(time_str.replace('Z', ''))
    else:
        dt = datetime.strptime(time_str, '%Y-%m-%d %H:%M:%S')
    return dt.year, dt.month, dt.day, dt.hour, dt.minute, dt.second


def classify_transient(ra, dec, anti_ra, anti_dec):
    """Classify a single transient as shadow/penumbra/sunlit."""
    sep = angular_separation(ra, dec, anti_ra, anti_dec)
    if sep <= GEO_SHADOW_RADIUS_DEG:
        return 'shadow'
    elif sep <= GEO_SHADOW_RADIUS_DEG + 1.0:
        return 'penumbra'
    else:
        return 'sunlit'


def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    repo_root = os.path.dirname(script_dir)
    vasco_path = os.path.join(repo_root, 'data', 'raw', 'vasco',
                              'SUPERVIKTIG_HELAVASCO_within2deg_CENTER.csv')
    if not os.path.exists(vasco_path):
        print(f"Data file not found: {vasco_path}")
        return

    plates = load_plates(vasco_path)
    n_plates = len(plates)
    total_transients = sum(len(m) for m in plates.values())
    print(f"Loaded {total_transients} transients across {n_plates} plates")
    print(f"Dataset: VASCO 2-degree center-of-plate filtered")
    print(f"Shadow cone radius: {GEO_SHADOW_RADIUS_DEG:.4f} deg")
    print(f"Plate field radius: {PLATE_FIELD_RADIUS_DEG:.1f} deg")

    total_overlap_area = 0.0
    total_plate_area = 0.0
    plate_area = math.pi * PLATE_FIELD_RADIUS_DEG ** 2  # sq deg per plate

    overlap_details = []
    shadow_transients = []
    penumbra_transients = []
    total_shadow = 0
    total_penumbra = 0
    total_sunlit = 0

    for plate_name, members in sorted(plates.items()):
        # 1. Plate center
        center_ra, center_dec = compute_plate_center(members)

        # 2. Observation time
        year, month, day, hour, minute, second = parse_obs_time(members[0]['time_str'])

        # 3. Anti-sun position
        anti_ra, anti_dec = anti_sun_position(year, month, day, hour, minute, second)

        # 4. Angular distance from plate center to shadow center
        sep = angular_separation(center_ra, center_dec, anti_ra, anti_dec)

        # 5. Overlap fraction
        frac = circle_overlap_fraction(sep, PLATE_FIELD_RADIUS_DEG, GEO_SHADOW_RADIUS_DEG)
        overlap = frac * plate_area

        total_overlap_area += overlap
        total_plate_area += plate_area

        # 6. Classify each transient individually
        plate_shadow = 0
        plate_penumbra = 0
        plate_sunlit = 0
        for t in members:
            cat = classify_transient(t['ra'], t['dec'], anti_ra, anti_dec)
            if cat == 'shadow':
                plate_shadow += 1
                shadow_transients.append({
                    'ra': t['ra'], 'dec': t['dec'], 'plate': plate_name,
                    'time': t['time_str'],
                    'sep_to_antisun': angular_separation(t['ra'], t['dec'], anti_ra, anti_dec)
                })
            elif cat == 'penumbra':
                plate_penumbra += 1
                penumbra_transients.append({
                    'ra': t['ra'], 'dec': t['dec'], 'plate': plate_name,
                    'time': t['time_str'],
                    'sep_to_antisun': angular_separation(t['ra'], t['dec'], anti_ra, anti_dec)
                })
            else:
                plate_sunlit += 1

        total_shadow += plate_shadow
        total_penumbra += plate_penumbra
        total_sunlit += plate_sunlit

        if frac > 0:
            overlap_details.append({
                'plate': plate_name,
                'center_ra': center_ra,
                'center_dec': center_dec,
                'anti_sun_ra': anti_ra,
                'anti_sun_dec': anti_dec,
                'separation': sep,
                'overlap_fraction': frac,
                'n_transients': len(members),
                'n_shadow': plate_shadow,
                'n_penumbra': plate_penumbra,
            })

    expected_shadow_fraction = total_overlap_area / total_plate_area
    observed_shadow_fraction = total_shadow / total_transients if total_transients > 0 else 0
    expected_shadow_transients = expected_shadow_fraction * total_transients

    print(f"\n{'='*70}")
    print(f"PLATE-BY-PLATE SHADOW COVERAGE RESULTS (2-deg Center-of-Plate Data)")
    print(f"{'='*70}")
    print(f"\n--- Dataset Summary ---")
    print(f"Total transients: {total_transients}")
    print(f"Number of plates: {n_plates}")
    print(f"Total plate area: {total_plate_area:.2f} sq deg")
    print(f"Total shadow overlap area: {total_overlap_area:.4f} sq deg")

    print(f"\n--- Expected vs Observed ---")
    print(f"Expected shadow fraction: {expected_shadow_fraction:.6f} "
          f"({100 * expected_shadow_fraction:.4f}%)")
    print(f"Expected shadow transients: {expected_shadow_transients:.1f}")
    print(f"Observed shadow transients: {total_shadow}")
    print(f"Observed shadow fraction: {observed_shadow_fraction:.6f} "
          f"({100 * observed_shadow_fraction:.4f}%)")

    if expected_shadow_transients > 0:
        deficit_ratio = expected_shadow_transients / max(total_shadow, 1)
        print(f"Deficit ratio: {deficit_ratio:.2f}x below expectation")

    print(f"\n--- Transient Classification ---")
    print(f"Shadow (within {GEO_SHADOW_RADIUS_DEG:.2f} deg of anti-sun): {total_shadow}")
    print(f"Penumbra (within +1 deg): {total_penumbra}")
    print(f"Sunlit: {total_sunlit}")

    print(f"\n--- Plates with Shadow Overlap ---")
    print(f"Plates overlapping shadow cone: {len(overlap_details)}")
    if overlap_details:
        n_full = sum(1 for d in overlap_details if d['overlap_fraction'] > 0.99)
        n_partial = len(overlap_details) - n_full
        print(f"  Fully inside shadow: {n_full}")
        print(f"  Partially overlapping: {n_partial}")

        print(f"\n{'Plate':<10} {'Sep(deg)':>10} {'Overlap%':>10} {'N_trans':>8} "
              f"{'N_shadow':>9} {'N_penumbra':>11}")
        print(f"{'-'*10} {'-'*10} {'-'*10} {'-'*8} {'-'*9} {'-'*11}")
        for d in sorted(overlap_details, key=lambda x: x['separation']):
            print(f"{d['plate']:<10} {d['separation']:>10.2f} "
                  f"{100*d['overlap_fraction']:>9.2f}% {d['n_transients']:>8} "
                  f"{d['n_shadow']:>9} {d['n_penumbra']:>11}")

    if shadow_transients:
        print(f"\n--- Individual Shadow Transients ---")
        print(f"{'RA':>10} {'Dec':>10} {'Plate':<10} {'Sep_antisun':>12} {'Time'}")
        for t in sorted(shadow_transients, key=lambda x: x['sep_to_antisun']):
            print(f"{t['ra']:>10.4f} {t['dec']:>10.4f} {t['plate']:<10} "
                  f"{t['sep_to_antisun']:>11.2f}° {t['time']}")

    # Save results
    out_dir = os.path.join(script_dir, 'results')
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, 'plate_shadow_coverage_2deg_center.txt')
    with open(out_path, 'w') as f:
        f.write(f"Plate-by-plate shadow coverage (2-deg center-of-plate data)\n")
        f.write(f"{'='*70}\n")
        f.write(f"Dataset: SUPERVIKTIG_HELAVASCO_within2deg_CENTER.csv\n")
        f.write(f"Shadow cone radius: {GEO_SHADOW_RADIUS_DEG:.4f} deg\n")
        f.write(f"Plate field radius: {PLATE_FIELD_RADIUS_DEG:.1f} deg\n")
        f.write(f"Total transients: {total_transients}\n")
        f.write(f"Number of plates: {n_plates}\n")
        f.write(f"Total plate area: {total_plate_area:.2f} sq deg\n")
        f.write(f"Total shadow overlap area: {total_overlap_area:.4f} sq deg\n\n")
        f.write(f"Expected shadow fraction: {expected_shadow_fraction:.6f} "
                f"({100 * expected_shadow_fraction:.4f}%)\n")
        f.write(f"Expected shadow transients: {expected_shadow_transients:.1f}\n")
        f.write(f"Observed shadow transients: {total_shadow}\n")
        f.write(f"Observed shadow fraction: {observed_shadow_fraction:.6f} "
                f"({100 * observed_shadow_fraction:.4f}%)\n")
        if expected_shadow_transients > 0:
            f.write(f"Deficit ratio: {expected_shadow_transients / max(total_shadow, 1):.2f}x\n")
        f.write(f"\nPlates with shadow overlap: {len(overlap_details)}\n")
        if overlap_details:
            for d in sorted(overlap_details, key=lambda x: x['separation']):
                f.write(f"  {d['plate']}: sep={d['separation']:.2f}°, "
                        f"overlap={100*d['overlap_fraction']:.2f}%, "
                        f"transients={d['n_transients']}, "
                        f"shadow={d['n_shadow']}, penumbra={d['n_penumbra']}\n")
        f.write(f"\nShadow transients: {total_shadow}\n")
        for t in sorted(shadow_transients, key=lambda x: x['sep_to_antisun']):
            f.write(f"  RA={t['ra']:.4f}, Dec={t['dec']:.4f}, plate={t['plate']}, "
                    f"sep={t['sep_to_antisun']:.2f}°, time={t['time']}\n")
    print(f"\nResults saved to {out_path}")

    return expected_shadow_fraction, observed_shadow_fraction


if __name__ == '__main__':
    main()
