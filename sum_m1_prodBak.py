#!/usr/bin/env python3

import numpy as np
import astropy.units as u
from astropy.table import Table
import warnings
from astropy.units import UnitsWarning


def summarize_ms1_product_sizes(tbl, cube_cap_GB=None):
    """
    Compute delivered product sizes for 'initial WSU' (MS1), then aggregate
    by project and by cycle.

    We interpret columns as follows (per row of the input table):
      - wsu_cubesize_initial_stepped2 [Gbyte]:
            size of ONE science target cube in initial WSU/MS1.
      - wsu_nspw_initial (dimensionless):
            number of SPWs we plan to deliver per target for MS1.
      - ntarget (dimensionless):
            number of targets in that project (doesn't change with WSU dev).
      - proposal_id [str]:
            ALMA project code.
      - cycle_info [str]:
            ALMA Cycle identifier.

    For each row:
      uncapped_row_prod  = cube_size * nspw * ntarget
      capped_row_prod    = min(cube_size, cube_cap_GB) * nspw * ntarget
        (if cube_cap_GB is None, no cap is applied)

    Then we:
      1. Sum rows with the same (proposal_id, cycle_info) to get project totals.
      2. Sum projects with the same cycle_info to get per-cycle totals.

    Returns
    -------
    proj_summary : list of dict
        Each dict has:
          {
            'proposal_id': str,
            'cycle_info': str,
            'prodsize_GB_capped': float,
            'prodsize_GB_uncapped': float,
          }

    cycle_summary : list of dict
        Each dict has:
          {
            'cycle_info': str,
            'total_prodsize_GB_capped': float,
            'total_prodsize_GB_uncapped': float,
          }
    """

    # Pull needed columns. We explicitly convert to numpy arrays / floats.
    cube_size = tbl['wsu_cubesize_initial_stepped2'].to(u.Gbyte)  # Quantity
    nspw      = np.asarray(tbl['wsu_nspw_initial'], dtype=float)
    ntarget   = np.asarray(tbl['ntarget'], dtype=float)

    proposal  = np.asarray(tbl['proposal_id']).astype(str)
    cycle     = np.asarray(tbl['cycle_info']).astype(str)

    # Uncapped per-row delivered product (Quantity in Gbyte)
    prod_uncap = cube_size * nspw * ntarget

    # Capped per-row delivered product (Quantity in Gbyte)
    if cube_cap_GB is not None:
        cap = cube_cap_GB * u.Gbyte
        eff_cube_size = np.minimum(cube_size, cap)  # elementwise Quantity min
    else:
        eff_cube_size = cube_size

    prod_cap = eff_cube_size * nspw * ntarget

    # -------------------------------------------------
    # 1. Accumulate per PROJECT (proposal_id, cycle_info)
    # -------------------------------------------------
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
        # Add (in Gbyte as plain floats)
        per_project[key]['uncap'] += unc_q.value
        per_project[key]['cap']   += cap_q.value

    # convert dict -> list of dicts for nicer downstream handling
    proj_summary = []
    for (pid, cyc), vals in per_project.items():
        proj_summary.append({
            'proposal_id': pid,
            'cycle_info': cyc,
            'prodsize_GB_capped': vals['cap'],
            'prodsize_GB_uncapped': vals['uncap'],
        })

    # -------------------------------------------------
    # 2. Accumulate per CYCLE across projects
    # -------------------------------------------------
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

    # sort cycle_summary by cycle_info just for pretty printing
    cycle_summary.sort(key=lambda d: d['cycle_info'])

    return proj_summary, cycle_summary


def main():
    # Astropy will complain about unknown units in unrelated columns
    # like 'Gvis' and 'Gvis / h'. We don't actually use those columns,
    # so let's just suppress that warning noise.
    warnings.simplefilter('ignore', UnitsWarning)

    # -------------------------------------------------
    # Load your table
    # NOTE the filename here matches what you showed via grep
    # -------------------------------------------------
    t = Table.read(
        "data/wsu_datarates_mit_per_mous_initial_goal_20250423.ecsv",
        format="ascii.ecsv",
    )

    # Case A: impose a 100 GB per-cube cap (the scenario you care about)
    proj_100, cycle_100 = summarize_ms1_product_sizes(t, cube_cap_GB=100.0)

    # Case B: NO cap at all (sanity check)
    proj_nocap, cycle_nocap = summarize_ms1_product_sizes(t, cube_cap_GB=None)

    # -------------------------------------------------
    # Print cycle-level comparison
    # -------------------------------------------------
    print("Per-cycle totals WITH 100 GB per-cube cap:")
    for row in cycle_100:
        print(
            f"{row['cycle_info']:>12}  "
            f"{row['total_prodsize_GB_capped']:12.1f} Gbyte"
        )

    print("\nPer-cycle totals with NO cap:")
    for row in cycle_nocap:
        print(
            f"{row['cycle_info']:>12}  "
            f"{row['total_prodsize_GB_uncapped']:12.1f} Gbyte"
        )

    # -------------------------------------------------
    # (Optional) If you want to see which projects are pigs:
    # sort projects by capped size desc and print top few.
    # -------------------------------------------------
    proj_sorted = sorted(
        proj_100,
        key=lambda d: d['prodsize_GB_capped'],
        reverse=True
    )
    print("\nTop 10 heaviest projects WITH 100 GB cap (GB):")
    for p in proj_sorted[:10]:
        print(
            f"{p['proposal_id']:>14}  {p['cycle_info']:>12}  "
            f"{p['prodsize_GB_capped']:12.1f} GB"
        )


if __name__ == "__main__":
    main()

