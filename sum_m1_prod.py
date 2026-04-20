#!/usr/bin/env python3

import numpy as np
import astropy.units as u
from astropy.table import Table
import warnings
from astropy.units import UnitsWarning

# -------------------------------------------------
# CONFIG
# -------------------------------------------------

# How many "cubes" per delivered SPW/target should we count?
# e.g. 2.0 = science cube + PB cube (or similar paired product)
#      1.0 = just the science cube
CUBE_SIZE_FACTOR = 2.0

# For human-readable printing at the end.
# Change to 1000.0 if you want decimal TB instead of TiB.
GB_PER_TB = 1024.0


def gb_to_tb(val_gb, gb_per_tb=GB_PER_TB):
    """Helper for nicer printing."""
    return val_gb / gb_per_tb


def summarize_ms1_product_sizes(tbl, cube_cap_GB=None, cube_size_factor=CUBE_SIZE_FACTOR):
    """
    Compute delivered product sizes for 'initial WSU' (MS1), aggregate by MOUS UID,
    then aggregate by cycle.

    Per input row we expect these columns:
      - wsu_cubesize_initial_stepped2 [Gbyte]:
            size of ONE science cube for this row's setup (MS1 assumptions).
      - wsu_nspw_initial:
            number of SPWs (i.e. distinct cubes) per target for MS1.
      - ntarget:
            number of targets for that MOUS/proposal.
      - mous:
            MOUS UID string.
      - cycle_info:
            ALMA Cycle label.

    The per-row product calculation is:

        uncapped_row_prod_GB
            = (cube_size_GB * cube_size_factor) * nspw * ntarget

        capped_row_prod_GB
            = (min(cube_size_GB, cube_cap_GB) * cube_size_factor)
              * nspw * ntarget

    where cube_size_factor (default 2.0) lets us account for
    "extra products per SPW" like PB cubes alongside the science cube.

    We then:
      1. Sum rows with the same (mous, cycle_info) to get one number per MOUS.
      2. Sum those MOUS totals per cycle_info to get per-cycle totals.

    Parameters
    ----------
    tbl : astropy.table.Table
    cube_cap_GB : float or None
        If not None, we clamp the *individual cube size* at this many GB
        BEFORE multiplying by nspw, ntarget, etc.
        This is the per-cube hard cap you asked for.
    cube_size_factor : float
        Multiplies the cube size before scaling up by SPWs and targets.
        Use 2.0 to include PB cubes, etc.

    Returns
    -------
    mous_summary : list of dict
        [
          {
            'mous': str,
            'cycle_info': str,
            'prodsize_GB_capped': float,
            'prodsize_GB_uncapped': float,
          },
          ...
        ]

    cycle_summary : list of dict
        [
          {
            'cycle_info': str,
            'total_prodsize_GB_capped': float,
            'total_prodsize_GB_uncapped': float,
          },
          ...
        ]
    """

    # pull needed columns
    cube_size = tbl['wsu_cubesize_initial_stepped2'].to(u.Gbyte)  # Quantity
    nspw      = np.asarray(tbl['wsu_nspw_initial'], dtype=float)
    ntarget   = np.asarray(tbl['ntarget'], dtype=float)

    mous_uid  = np.asarray(tbl['mous']).astype(str)
    cycle     = np.asarray(tbl['cycle_info']).astype(str)

    # per-row UNCAPPED product (Quantity in Gbyte)
    prod_uncap = cube_size * cube_size_factor * nspw * ntarget

    # per-row CAPPED product (Quantity in Gbyte)
    if cube_cap_GB is not None:
        cap = cube_cap_GB * u.Gbyte
        eff_cube_size = np.minimum(cube_size, cap)  # Quantity min
    else:
        eff_cube_size = cube_size

    prod_cap = eff_cube_size * cube_size_factor * nspw * ntarget

    # -------------------------------------------------
    # 1. accumulate per (mous, cycle_info)
    # -------------------------------------------------
    per_mous = {}
    for mu, cyc, unc_q, cap_q in zip(
        mous_uid,
        cycle,
        prod_uncap.to(u.Gbyte),
        prod_cap.to(u.Gbyte),
    ):
        key = (mu, cyc)
        if key not in per_mous:
            per_mous[key] = {
                'uncap': 0.0,
                'cap':   0.0,
            }
        per_mous[key]['uncap'] += unc_q.value  # float GB
        per_mous[key]['cap']   += cap_q.value  # float GB

    mous_summary = []
    for (mu, cyc), vals in per_mous.items():
        mous_summary.append({
            'mous': mu,
            'cycle_info': cyc,
            'prodsize_GB_capped': vals['cap'],
            'prodsize_GB_uncapped': vals['uncap'],
        })

    # -------------------------------------------------
    # 2. accumulate per cycle from the MOUS-level totals
    # -------------------------------------------------
    per_cycle = {}
    for m in mous_summary:
        cyc = m['cycle_info']
        if cyc not in per_cycle:
            per_cycle[cyc] = {
                'uncap': 0.0,
                'cap':   0.0,
            }
        per_cycle[cyc]['uncap'] += m['prodsize_GB_uncapped']
        per_cycle[cyc]['cap']   += m['prodsize_GB_capped']

    cycle_summary = []
    for cyc, vals in per_cycle.items():
        cycle_summary.append({
            'cycle_info': cyc,
            'total_prodsize_GB_capped': vals['cap'],
            'total_prodsize_GB_uncapped': vals['uncap'],
        })

    # stable order by cycle label for printing
    cycle_summary.sort(key=lambda d: d['cycle_info'])

    return mous_summary, cycle_summary


def build_mous_details(mous_list, tbl):
    """
    Enrich each MOUS entry with:
      - representative cube size in M1 (GB),
      - representative wsu_nspw_initial,
      - representative ntarget,
      - capped and uncapped totals (GB),
      - plus TB conversions.

    We do this by re-querying the original table rows for that (mous, cycle_info).
    """
    details = []

    for m in mous_list:
        mu = m['mous']
        cyc = m['cycle_info']

        mask = (
            (tbl['mous'].astype(str) == mu) &
            (tbl['cycle_info'].astype(str) == cyc)
        )
        rows = tbl[mask]

        if len(rows) == 0:
            cube_gb = np.nan
            nspw_proj = np.nan
            ntargets_proj = np.nan
        else:
            # Representative stats per MOUS:
            # - cube size: mean of wsu_cubesize_initial_stepped2 (in GB)
            # - NSPW and Ntargets: max across matching rows (robust to repeats)
            cube_gb = np.mean(
                rows['wsu_cubesize_initial_stepped2'].to(u.Gbyte)
            ).value

            nspw_proj = float(
                np.max(np.asarray(rows['wsu_nspw_initial'], dtype=float))
            )

            ntargets_proj = float(
                np.max(np.asarray(rows['ntarget'], dtype=float))
            )

        capped_gb   = m['prodsize_GB_capped']
        uncapped_gb = m['prodsize_GB_uncapped']

        details.append({
            'mous': mu,
            'cycle_info': cyc,
            'cube_gb': cube_gb,
            'nspw': nspw_proj,
            'ntarget': ntargets_proj,
            'prod_cap_gb': capped_gb,
            'prod_uncap_gb': uncapped_gb,
            'prod_cap_tb': gb_to_tb(capped_gb),
            'prod_uncap_tb': gb_to_tb(uncapped_gb),
        })

    return details


def print_top10(title, details, sort_key):
    """
    Print a Top 10 list with aligned columns.
    sort_key: a lambda that takes a detail dict and returns the sort value.
    """
    print(f"\n*** {title} ***")

    # sort descending by the provided key
    sorted_details = sorted(details, key=sort_key, reverse=True)

    # consistent column widths between header and rows
    # 30s  8s  14s 14s 14s 7s 10s
    header = (
        f"{'MOUSuid':>30}  "
        f"{'Cycle':>8}  "
        f"{'ProdCap[TB]':>14}  "
        f"{'ProdUncap[TB]':>14}  "
        f"{'CubeM1[GB]':>14}  "
        f"{'NSPW':>7}  "
        f"{'Ntargets':>10}"
    )
    print(header)

    for d in sorted_details[:10]:
        print(
            f"{d['mous']:>30}  "
            f"{d['cycle_info']:>8}  "
            f"{d['prod_cap_tb']:14.2f}  "
            f"{d['prod_uncap_tb']:14.2f}  "
            f"{d['cube_gb']:14.2f}  "
            f"{d['nspw']:7.0f}  "
            f"{d['ntarget']:10.0f}"
        )


def main():
    # suppress warnings about weird units ('Gvis', etc.)
    warnings.simplefilter('ignore', UnitsWarning)

    # -------------------------------------------------
    # Load table
    # NOTE: this filename matches what you've been using
    # -------------------------------------------------
    t = Table.read(
        "data/wsu_datarates_mit_per_mous_initial_goal_20250423.ecsv",
        format="ascii.ecsv",
    )

    # with 100 GB per-cube cap
    mous_100, cycle_100 = summarize_ms1_product_sizes(
        t,
        cube_cap_GB=100.0,
        cube_size_factor=CUBE_SIZE_FACTOR,
    )

    # no cap (sanity check)
    mous_nocap, cycle_nocap = summarize_ms1_product_sizes(
        t,
        cube_cap_GB=None,
        cube_size_factor=CUBE_SIZE_FACTOR,
    )

    # -------------------------------------------------
    # Print cycle-level totals in TB
    # ALSO include uncapped totals now
    # -------------------------------------------------
    print(f"Per-cycle totals WITH 100 GB per-cube cap (TB), "
          f"CUBE_SIZE_FACTOR={CUBE_SIZE_FACTOR}:")
    for row in cycle_100:
        tb_cap   = gb_to_tb(row['total_prodsize_GB_capped'])
        tb_uncap = gb_to_tb(row['total_prodsize_GB_uncapped'])
        print(
            f"{row['cycle_info']:>12}  "
            f"{tb_cap:12.2f} Tbyte (capped)   "
            f"{tb_uncap:12.2f} Tbyte (uncapped)"
        )

    print(f"\nPer-cycle totals with NO cap (TB), "
          f"CUBE_SIZE_FACTOR={CUBE_SIZE_FACTOR}:")
    for row in cycle_nocap:
        tb_uncap = gb_to_tb(row['total_prodsize_GB_uncapped'])
        print(
            f"{row['cycle_info']:>12}  "
            f"{tb_uncap:12.2f} Tbyte (uncapped)"
        )

    # -------------------------------------------------
    # Build enriched per-MOUS details (capped run),
    # because we want to sort/print 3 different ways
    # -------------------------------------------------
    mous_details = build_mous_details(mous_100, t)

    # 1. Top 10 by UN-CAPPED aggregate product size
    print_top10(
        "Top 10 MOUS by UN-CAPPED aggregate product size (largest total TB)",
        mous_details,
        sort_key=lambda d: d['prod_uncap_tb'],
    )

    # 2. Top 10 by cube size
    print_top10(
        "Top 10 MOUS by cube size (largest individual cube in GB, averaged across rows)",
        mous_details,
        sort_key=lambda d: d['cube_gb'],
    )

    # 3. Top 10 by CAPPED aggregate product size
    print_top10(
        "Top 10 MOUS by CAPPED aggregate product size (largest total TB after per-cube cap)",
        mous_details,
        sort_key=lambda d: d['prod_cap_tb'],
    )


if __name__ == "__main__":
    main()

