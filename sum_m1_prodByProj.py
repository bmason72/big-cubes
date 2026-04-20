#!/usr/bin/env python3

import numpy as np
import astropy.units as u
from astropy.table import Table
import warnings
from astropy.units import UnitsWarning


GB_PER_TB = 1024.0  # change to 1000.0 if you want decimal TB


def summarize_ms1_product_sizes(tbl, cube_cap_GB=None):
    """
    Compute delivered product sizes for 'initial WSU' (MS1), then aggregate
    by project and by cycle.

    Per-row interpretation:
      - wsu_cubesize_initial_stepped2 [Gbyte]:
            size of ONE science-target, single-SPW cube for MS1.
      - wsu_nspw_initial:
            how many SPW cubes we deliver per target.
      - ntarget:
            number of targets in that project.
      - proposal_id:
            ALMA project code.
      - cycle_info:
            ALMA Cycle label.

    Row math:
        uncapped_row_prod_GB = cube_size_GB * nspw * ntarget
        capped_row_prod_GB   = min(cube_size_GB, cube_cap_GB) * nspw * ntarget
                               (if cube_cap_GB is None, no cap)

    Then:
      - sum rows with the same (proposal_id, cycle_info) → project totals
      - sum project totals with the same cycle_info     → cycle totals

    Returns
    -------
    proj_summary : list of dict
        [
          {
            'proposal_id': str,
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

    # pull columns
    cube_size = tbl['wsu_cubesize_initial_stepped2'].to(u.Gbyte)  # Quantity
    nspw      = np.asarray(tbl['wsu_nspw_initial'], dtype=float)
    ntarget   = np.asarray(tbl['ntarget'], dtype=float)

    proposal  = np.asarray(tbl['proposal_id']).astype(str)
    cycle     = np.asarray(tbl['cycle_info']).astype(str)

    # per-row uncapped product (Quantity in Gbyte)
    prod_uncap = cube_size * nspw * ntarget

    # per-row capped product (Quantity in Gbyte)
    if cube_cap_GB is not None:
        cap = cube_cap_GB * u.Gbyte
        eff_cube_size = np.minimum(cube_size, cap)  # elementwise Quantity min
    else:
        eff_cube_size = cube_size

    prod_cap = eff_cube_size * nspw * ntarget

    # accumulate per (proposal_id, cycle_info)
    per_project = {}
    for pid, cyc, unc_q, cap_q in zip(
        proposal,
        cycle,
        prod_uncap.to(u.Gbyte),
        prod_cap.to(u.Gbyte),
    ):
        key = (pid, cyc)
        if key not in per_project:
            per_project[key] = {
                'uncap': 0.0,
                'cap':   0.0,
            }
        per_project[key]['uncap'] += unc_q.value  # float GB
        per_project[key]['cap']   += cap_q.value  # float GB

    proj_summary = []
    for (pid, cyc), vals in per_project.items():
        proj_summary.append({
            'proposal_id': pid,
            'cycle_info': cyc,
            'prodsize_GB_capped': vals['cap'],
            'prodsize_GB_uncapped': vals['uncap'],
        })

    # accumulate per cycle
    per_cycle = {}
    for proj in proj_summary:
        cyc = proj['cycle_info']
        if cyc not in per_cycle:
            per_cycle[cyc] = {
                'uncap': 0.0,
                'cap':   0.0,
            }
        per_cycle[cyc]['uncap'] += proj['prodsize_GB_uncapped']
        per_cycle[cyc]['cap']   += proj['prodsize_GB_capped']

    cycle_summary = []
    for cyc, vals in per_cycle.items():
        cycle_summary.append({
            'cycle_info': cyc,
            'total_prodsize_GB_capped': vals['cap'],
            'total_prodsize_GB_uncapped': vals['uncap'],
        })

    # sort cycle output for stable printing
    cycle_summary.sort(key=lambda d: d['cycle_info'])

    return proj_summary, cycle_summary


def main():
    # suppress unknown-unit spam like 'Gvis'
    warnings.simplefilter('ignore', UnitsWarning)

    # load the table from your working path
    t = Table.read(
        "data/wsu_datarates_mit_per_mous_initial_goal_20250423.ecsv",
        format="ascii.ecsv",
    )

    # with 100 GB hard per-cube cap
    proj_100, cycle_100 = summarize_ms1_product_sizes(t, cube_cap_GB=100.0)

    # no cap (sanity check)
    proj_nocap, cycle_nocap = summarize_ms1_product_sizes(t, cube_cap_GB=None)

    # helper for GB → TB
    def gb_to_tb(val_gb):
        return val_gb / GB_PER_TB

    # print per-cycle totals
    print("Per-cycle totals WITH 100 GB per-cube cap (TB):")
    for row in cycle_100:
        tb_cap = gb_to_tb(row['total_prodsize_GB_capped'])
        print(
            f"{row['cycle_info']:>12}  "
            f"{tb_cap:12.2f} Tbyte"
        )

    print("\nPer-cycle totals with NO cap (TB):")
    for row in cycle_nocap:
        tb_uncap = gb_to_tb(row['total_prodsize_GB_uncapped'])
        print(
            f"{row['cycle_info']:>12}  "
            f"{tb_uncap:12.2f} Tbyte"
        )

    # biggest projects by capped size
    proj_sorted = sorted(
        proj_100,
        key=lambda d: d['prodsize_GB_capped'],
        reverse=True
    )

    print("\nTop 10 heaviest projects WITH 100 GB cap (TB):")
    for p in proj_sorted[:10]:
        tb_cap_proj = gb_to_tb(p['prodsize_GB_capped'])
        print(
            f"{p['proposal_id']:>14}  {p['cycle_info']:>12}  "
            f"{tb_cap_proj:12.2f} TB"
        )


if __name__ == "__main__":
    main()

