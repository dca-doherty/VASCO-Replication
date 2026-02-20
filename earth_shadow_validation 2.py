#!/usr/bin/env python3
"""
Earth Shadow Deficit Validation
Replicates methodology from Villarroel et al. PASP paper.

Objects detected inside Earth's geometric shadow cannot be sunlit debris
since no solar illumination reaches the umbra region. This provides a
geometric test for the transient population.

Author: Brian Doherty
Date: December 2025
Contact: briandohertyresearch@gmail.com
"""

import csv
import math
from datetime import datetime
import os

# physical constants
EARTH_RADIUS_KM = 6371.0
SUN_RADIUS_KM = 696000.0
AU_KM = 149597870.7

# GEO orbital parameters
GEO_ORBITAL_RADIUS_KM = 42164.0

# umbra geometry calculations
# umbra length is distance from Earth center to cone vertex
UMBRA_LENGTH_KM = EARTH_RADIUS_KM * AU_KM / (SUN_RADIUS_KM - EARTH_RADIUS_KM)

# at GEO distance, umbra radius is reduced due to cone narrowing
UMBRA_RADIUS_AT_GEO = EARTH_RADIUS_KM * (1 - GEO_ORBITAL_RADIUS_KM / UMBRA_LENGTH_KM)

# shadow angle subtended at GEO
GEO_SHADOW_RADIUS_DEG = math.degrees(math.asin(UMBRA_RADIUS_AT_GEO / GEO_ORBITAL_RADIUS_KM))

# simple geometric approximation for comparison
SIMPLE_SHADOW_RADIUS_DEG = math.degrees(math.asin(EARTH_RADIUS_KM / GEO_ORBITAL_RADIUS_KM))


def sun_position_meeus(year, month, day, hour=12, minute=0, second=0):
    """
    Calculate Sun RA and Dec using Meeus algorithm.
    Based on Astronomical Algorithms by Jean Meeus (2nd edition).
    Accuracy is approximately 1 arcminute.
    """
    y, m = year, month
    if m <= 2:
        y -= 1
        m += 12

    # julian date
    A = int(y / 100)
    B = 2 - A + int(A / 4)
    JD = int(365.25 * (y + 4716)) + int(30.6001 * (m + 1)) + day + B - 1524.5
    JD += (hour + minute / 60 + second / 3600) / 24.0

    # julian centuries from J2000.0
    T = (JD - 2451545.0) / 36525.0

    # geometric mean longitude of sun
    L0 = (280.46646 + 36000.76983 * T + 0.0003032 * T ** 2) % 360

    # mean anomaly
    M = (357.52911 + 35999.05029 * T - 0.0001537 * T ** 2) % 360
    M_rad = math.radians(M)

    # eccentricity
    e = 0.016708634 - 0.000042037 * T - 0.0000001267 * T ** 2

    # equation of center
    C = ((1.914602 - 0.004817 * T - 0.000014 * T ** 2) * math.sin(M_rad)
         + (0.019993 - 0.000101 * T) * math.sin(2 * M_rad)
         + 0.000289 * math.sin(3 * M_rad))

    # true longitude
    sun_lon = L0 + C

    # apparent longitude with nutation and aberration
    omega = 125.04 - 1934.136 * T
    sun_lon_app = sun_lon - 0.00569 - 0.00478 * math.sin(math.radians(omega))
    sun_lon_rad = math.radians(sun_lon_app)

    # obliquity of ecliptic
    eps0 = 23.439291111 - 0.013004167 * T - 0.00000016389 * T ** 2 + 0.0000005036 * T ** 3
    eps = eps0 + 0.00256 * math.cos(math.radians(omega))
    eps_rad = math.radians(eps)

    # equatorial coordinates
    ra_rad = math.atan2(math.cos(eps_rad) * math.sin(sun_lon_rad), math.cos(sun_lon_rad))
    dec_rad = math.asin(math.sin(eps_rad) * math.sin(sun_lon_rad))

    return math.degrees(ra_rad) % 360, math.degrees(dec_rad)


def anti_sun_position(year, month, day, hour=12, minute=0, second=0):
    """Calculate anti-sun position (shadow center)."""
    sun_ra, sun_dec = sun_position_meeus(year, month, day, hour, minute, second)
    anti_ra = (sun_ra + 180.0) % 360
    anti_dec = -sun_dec
    return anti_ra, anti_dec


def angular_separation(ra1, dec1, ra2, dec2):
    """Calculate angular separation using haversine formula."""
    ra1, dec1, ra2, dec2 = map(math.radians, [ra1, dec1, ra2, dec2])
    dra = ra2 - ra1
    ddec = dec2 - dec1
    a = math.sin(ddec / 2) ** 2 + math.cos(dec1) * math.cos(dec2) * math.sin(dra / 2) ** 2
    return math.degrees(2 * math.asin(math.sqrt(min(1, max(0, a)))))


def load_vasco_transients(filepath):
    """Load VASCO transient catalog."""
    transients = []
    errors = 0

    with open(filepath, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                ra = float(row['RA'])
                dec = float(row['Dec'])
                name = row.get('Name', '')
                time_str = row['UTobservation']

                # parse datetime
                if 'T' in time_str:
                    dt = datetime.fromisoformat(time_str.replace('Z', ''))
                else:
                    dt = datetime.strptime(time_str, '%Y-%m-%d %H:%M:%S')

                transients.append({
                    'name': name,
                    'ra': ra,
                    'dec': dec,
                    'date': dt
                })

            except Exception:
                errors += 1

    return transients, errors


def classify_shadow_transients(transients):
    """Classify transients by shadow position."""
    results = []

    for t in transients:
        dt = t['date']
        anti_ra, anti_dec = anti_sun_position(
            dt.year, dt.month, dt.day,
            dt.hour, dt.minute, dt.second
        )

        shadow_dist = angular_separation(t['ra'], t['dec'], anti_ra, anti_dec)
        in_umbra = shadow_dist <= GEO_SHADOW_RADIUS_DEG
        in_simple_shadow = shadow_dist <= SIMPLE_SHADOW_RADIUS_DEG

        results.append({
            'name': t['name'],
            'ra': t['ra'],
            'dec': t['dec'],
            'date': t['date'],
            'shadow_distance': shadow_dist,
            'in_umbra': in_umbra,
            'in_simple_shadow': in_simple_shadow,
            'anti_sun_ra': anti_ra,
            'anti_sun_dec': anti_dec
        })

    return results


def run_shadow_analysis(vasco_path, output_dir=None):
    """Main shadow analysis routine."""
    if output_dir is None:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        output_dir = os.path.join(script_dir, 'results')
    os.makedirs(output_dir, exist_ok=True)

    print("Shadow geometry parameters:")
    print(f"  Earth radius: {EARTH_RADIUS_KM:,.0f} km")
    print(f"  GEO orbital radius: {GEO_ORBITAL_RADIUS_KM:,.0f} km")
    print(f"  Umbra length: {UMBRA_LENGTH_KM:,.0f} km")
    print(f"  Umbra radius at GEO: {UMBRA_RADIUS_AT_GEO:,.1f} km")
    print(f"  Shadow angle (proper): {GEO_SHADOW_RADIUS_DEG:.4f} deg")
    print(f"  Shadow angle (simple): {SIMPLE_SHADOW_RADIUS_DEG:.4f} deg")

    # load transients
    print(f"\nLoading VASCO data from: {vasco_path}")
    transients, errors = load_vasco_transients(vasco_path)
    print(f"  Loaded: {len(transients):,} transients")
    if errors > 0:
        print(f"  Parse errors: {errors}")

    # classify
    print("\nClassifying transients by shadow position...")
    results = classify_shadow_transients(transients)

    umbra_count = sum(1 for r in results if r['in_umbra'])
    simple_count = sum(1 for r in results if r['in_simple_shadow'])

    print(f"\nResults:")
    print(f"  Total transients: {len(results):,}")
    print(f"  In proper umbra (< {GEO_SHADOW_RADIUS_DEG:.2f} deg): {umbra_count:,} ({100 * umbra_count / len(results):.3f}%)")
    print(f"  In simple shadow (< {SIMPLE_SHADOW_RADIUS_DEG:.2f} deg): {simple_count:,} ({100 * simple_count / len(results):.3f}%)")

    # save results
    output_path = os.path.join(output_dir, 'shadow_classification.csv')
    with open(output_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Name', 'RA', 'Dec', 'Date', 'shadow_distance',
                         'in_umbra', 'in_simple_shadow'])
        for r in results:
            writer.writerow([
                r['name'], r['ra'], r['dec'],
                r['date'].strftime('%Y-%m-%d %H:%M:%S'),
                r['shadow_distance'], r['in_umbra'], r['in_simple_shadow']
            ])

    print(f"\nResults saved to {output_path}")

    # save umbra-only subset
    umbra_transients = [r for r in results if r['in_umbra']]
    umbra_path = os.path.join(output_dir, 'umbra_transients.csv')
    with open(umbra_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Name', 'RA', 'Dec', 'Date', 'shadow_distance'])
        for r in umbra_transients:
            writer.writerow([
                r['name'], r['ra'], r['dec'],
                r['date'].strftime('%Y-%m-%d %H:%M:%S'),
                r['shadow_distance']
            ])

    print(f"Umbra transients saved to {umbra_path}")

    return {
        'total': len(results),
        'umbra_count': umbra_count,
        'simple_count': simple_count,
        'umbra_fraction': umbra_count / len(results),
        'results': results
    }


if __name__ == '__main__':
    # look for data file in local data folder first
    script_dir = os.path.dirname(os.path.abspath(__file__))
    vasco_paths = [
        os.path.join(script_dir, 'data', 'SUPERVIKTIG_HELAVASCO.csv'),
        'SUPERVIKTIG_HELAVASCO.csv',
        'data/SUPERVIKTIG_HELAVASCO.csv'
    ]

    vasco_path = None
    for path in vasco_paths:
        if os.path.exists(path):
            vasco_path = path
            break

    if vasco_path:
        results = run_shadow_analysis(vasco_path)
    else:
        print("VASCO data file not found")
