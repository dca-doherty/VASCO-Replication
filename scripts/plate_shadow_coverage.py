#!/usr/bin/env python3
"""
Plate-by-plate shadow coverage calculation.

For each POSS-I plate, compute:
  1. Plate center (mean RA/Dec of transients on the plate)
  2. Observation time from plate metadata
  3. Anti-sun position at that time
  4. Fraction of the plate's 2-degree circular field overlapping the shadow cone
  5. Sum to get total expected shadow coverage as a fraction of total plate area

This replaces the naive full-sky geometric expectation (~1.4%) with a precise
plate-by-plate calculation that accounts for actual sky coverage.

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
    # Using the standard lens area formula
    cos_arg1 = (d * d + r1 * r1 - r2 * r2) / (2 * d * r1)
    cos_arg2 = (d * d + r2 * r2 - r1 * r1) / (2 * d * r2)
    # clamp for numerical safety
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


def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    vasco_path = os.path.join(script_dir, 'data', 'SUPERVIKTIG_HELAVASCO.csv')
    if not os.path.exists(vasco_path):
        print(f"Data file not found: {vasco_path}")
        return

    plates = load_plates(vasco_path)
    n_plates = len(plates)
    print(f"Loaded {n_plates} plates")
    print(f"Shadow cone radius: {GEO_SHADOW_RADIUS_DEG:.4f} deg")
    print(f"Plate field radius: {PLATE_FIELD_RADIUS_DEG:.1f} deg")

    total_overlap_area = 0.0
    total_plate_area = 0.0
    plate_area = math.pi * PLATE_FIELD_RADIUS_DEG ** 2  # sq deg per plate

    overlap_details = []

    for plate_name, members in sorted(plates.items()):
        # 1. Plate center
        center_ra, center_dec = compute_plate_center(members)

        # 2. Observation time (all transients on a plate share the same time)
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
            })

    expected_shadow_fraction = total_overlap_area / total_plate_area

    print(f"\n{'='*60}")
    print(f"PLATE-BY-PLATE SHADOW COVERAGE RESULTS")
    print(f"{'='*60}")
    print(f"Number of plates: {n_plates}")
    print(f"Total plate area: {total_plate_area:.2f} sq deg")
    print(f"Total shadow overlap area: {total_overlap_area:.4f} sq deg")
    print(f"Expected shadow fraction: {expected_shadow_fraction:.6f} "
          f"({100 * expected_shadow_fraction:.4f}%)")
    print(f"\nPlates with any shadow overlap: {len(overlap_details)}")

    if overlap_details:
        print(f"\nPlates overlapping shadow cone:")
        print(f"{'Plate':<10} {'Sep(deg)':>10} {'Overlap%':>10} {'N_trans':>8}")
        for d in sorted(overlap_details, key=lambda x: x['separation']):
            print(f"{d['plate']:<10} {d['separation']:>10.2f} "
                  f"{100*d['overlap_fraction']:>9.2f}% {d['n_transients']:>8}")

    # Save results
    out_dir = os.path.join(script_dir, 'results')
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, 'plate_shadow_coverage.txt')
    with open(out_path, 'w') as f:
        f.write(f"Plate-by-plate shadow coverage calculation\n")
        f.write(f"Shadow cone radius: {GEO_SHADOW_RADIUS_DEG:.4f} deg\n")
        f.write(f"Plate field radius: {PLATE_FIELD_RADIUS_DEG:.1f} deg\n")
        f.write(f"Number of plates: {n_plates}\n")
        f.write(f"Total plate area: {total_plate_area:.2f} sq deg\n")
        f.write(f"Total shadow overlap area: {total_overlap_area:.4f} sq deg\n")
        f.write(f"Expected shadow fraction: {expected_shadow_fraction:.6f} "
                f"({100 * expected_shadow_fraction:.4f}%)\n")
        f.write(f"Plates with shadow overlap: {len(overlap_details)}\n")
    print(f"\nResults saved to {out_path}")

    return expected_shadow_fraction


if __name__ == '__main__':
    main()
