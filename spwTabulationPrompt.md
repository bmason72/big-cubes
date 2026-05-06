Here is the full updated prompt.

```text
Write a Python script to parse a collection of ALMA Cycle 11 pipeline weblog `listobs` text files, accompanied where available by `cont.dat` files and an optional metadata CSV.

Goal:
Create a reliable visibility spectral-setup catalog for ALMA MOUS/SB/EB observations, suitable for summarizing and grouping PI spectral configuration choices.

The core output should support:
1. Estimating the native visibility spectral setup for each MOUS.
2. Determining NSPW distributions.
3. Measuring per-SPW channel count, bandwidth, channel width, and velocity resolution.
4. Building compact setup signatures based on binned SPW bandwidth and velocity resolution.
5. Using `cont.dat` to estimate what fraction of each SPW is associated with line-contaminated/non-continuum frequency ranges.
6. Producing both machine-readable tables and human-readable summary statistics.

Directory/file conventions:
1. The top-level directory name encodes the MOUS UID.
2. Each MOUS directory may contain N_EB separate `listobs` files.
3. The `listobs` filename encodes the EB UID.
4. Spectral setups in different EBs within the same MOUS should usually be very similar.
5. Parse all EBs for sanity checking, but allow MOUS-level aggregation using:
   - the first EB as the representative setup, if all EBs are consistent;
   - or an aggregated/consensus setup if small differences are detected.
6. Preserve provenance:
   - top-level MOUS directory name
   - parsed MOUS UID
   - listobs filename
   - parsed EB UID
   - cont.dat filename, if present
   - parse warnings/errors

Inputs:
1. Directory tree containing one top-level directory per MOUS.
2. Each MOUS directory contains one or more `listobs` files, one per EB.
3. Optional matching `cont.dat` files.
4. Optional metadata CSV with columns such as:
   - project_code
   - member_ous_uid
   - group_ous_uid
   - science_goal_uid
   - schedblock_name
   - science_category or scientific_category
   - science_keyword
   - band
   - array_type
   - config_proxy
   - observation_date
   - eb_uid / asdm_uid
   - listobs_path
   - contdat_path

For each `listobs` file:
1. Parse the spectral-window table, extracting:
   - spw_id
   - spw_name
   - nchan
   - frame
   - ch0 frequency
   - channel width
   - total bandwidth
   - correlation products
2. Parse the field/source table:
   - field_id
   - field_name / target_name
3. Parse the scan table:
   - scan_id
   - field_id / field_name
   - scan intent
   - spw_ids used in that scan
4. Identify science-target scans, primarily those with scan intent containing `OBSERVE_TARGET`.
5. For each science target, determine the science SPWs actually used by joining target scans to the SPW table.
6. Exclude obvious auxiliary/non-science SPWs such as WVR, calibration-only, or SPWs not used in science-target scans.
7. For each target/setup/SPW, compute:
   - center_freq_hz
   - freq_min_hz
   - freq_max_hz
   - total_bw_hz
   - abs_chan_width_hz
   - nchan
   - velocity_resolution_kms = c * abs_chan_width_hz / center_freq_hz
   - polarization/correlation products

For each `cont.dat` file:
1. Parse all target-specific and SPW-specific continuum ranges.
2. Treat the ranges explicitly listed in `cont.dat` as continuum ranges.
3. Treat the complement of the continuum ranges within each SPW as the line ranges.
4. For each target/SPW, compute the union of all non-continuum frequency intervals.
5. Compute two line-coverage metrics:

   A. Span-based line range:
      - line_range_min_hz = min(line-range frequency)
      - line_range_max_hz = max(line-range frequency)
      - line_range_width_hz = line_range_max_hz - line_range_min_hz

   B. Summed disjoint line coverage:
      - summed_line_width_hz = sum(width of each disjoint non-continuum interval)

6. Add a buffer of 10 channels on either side of each line metric:

   A. For the span-based line range:
      - buffered_line_range_min_hz = line_range_min_hz - 10 * abs_chan_width_hz
      - buffered_line_range_max_hz = line_range_max_hz + 10 * abs_chan_width_hz
      - clip/guard against overflowing the SPW band edges
      - buffered_line_range_width_hz = buffered_line_range_max_hz - buffered_line_range_min_hz

   B. For the disjoint line intervals:
      - expand each disjoint non-continuum interval by 10 channels on either side
      - clip/guard against overflowing the SPW band edges
      - merge overlapping buffered intervals
      - buffered_summed_line_width_hz = sum(width of each merged buffered interval)

7. For each target/SPW, compute:
   - line_range_min_hz
   - line_range_max_hz
   - buffered_line_range_min_hz
   - buffered_line_range_max_hz
   - line_range_width_hz
   - buffered_line_range_width_hz
   - summed_line_width_hz
   - buffered_summed_line_width_hz
   - line_range_fraction_of_spw = line_range_width_hz / total_bw_hz
   - buffered_line_range_fraction_of_spw = buffered_line_range_width_hz / total_bw_hz
   - summed_line_fraction_of_spw = summed_line_width_hz / total_bw_hz
   - buffered_summed_line_fraction_of_spw = buffered_summed_line_width_hz / total_bw_hz
8. `cont.dat` may contain multiple targets and possibly multiple EB-like sections. Preserve target/EB provenance if present.
9. For MOUS-level summaries, average line-range fractions over targets and over EBs if applicable.
10. Report both:
   - per-target/per-SPW line-range fractions
   - MOUS-level average line-range fraction per SPW
11. If a SPW has no apparent line range after taking the complement of continuum ranges, set the line-range fractions to 0 and flag it.
12. If a SPW has no usable cont.dat information, set values to NaN and flag it.
13. Store both span-based and summed-disjoint metrics. The span-based metric captures the total frequency extent over which line contamination occurs, while the summed-disjoint metric captures the actual fractional coverage by non-continuum intervals.

Outputs:

1. `spw_level_table.csv`

   One row per science-target SPW per EB, with visibility spectral metadata:
   - mous_uid
   - eb_uid
   - listobs_path
   - target_name
   - field_id
   - spw_id
   - spw_name
   - nchan
   - center_freq_hz
   - freq_min_hz
   - freq_max_hz
   - total_bw_hz
   - abs_chan_width_hz
   - velocity_resolution_kms
   - corr_products
   - science_target_flag
   - parse_warning

2. `contdat_line_ranges.csv`

   One row per target/SPW/cont.dat entry:
   - mous_uid
   - eb_uid, if available
   - contdat_path
   - target_name
   - spw_id
   - continuum_ranges_hz
   - disjoint_line_ranges_hz
   - buffered_disjoint_line_ranges_hz
   - line_range_min_hz
   - line_range_max_hz
   - buffered_line_range_min_hz
   - buffered_line_range_max_hz
   - line_range_width_hz
   - buffered_line_range_width_hz
   - summed_line_width_hz
   - buffered_summed_line_width_hz
   - line_range_fraction_of_spw
   - buffered_line_range_fraction_of_spw
   - summed_line_fraction_of_spw
   - buffered_summed_line_fraction_of_spw
   - contdat_parse_warning

3. `mous_spw_summary.csv`

   One row per MOUS/SPW after aggregating over EBs and targets:
   - mous_uid
   - representative_eb_uid
   - spw_id
   - nchan
   - center_freq_hz
   - freq_min_hz
   - freq_max_hz
   - total_bw_hz
   - abs_chan_width_hz
   - velocity_resolution_kms
   - mean_line_range_fraction_of_spw
   - mean_buffered_line_range_fraction_of_spw
   - mean_summed_line_fraction_of_spw
   - mean_buffered_summed_line_fraction_of_spw
   - n_targets_with_contdat
   - n_ebs_parsed
   - eb_consistency_flag
   - eb_consistency_notes

4. `setup_level_table.csv`

   One row per MOUS or MOUS/target setup:
   - mous_uid
   - representative_eb_uid
   - n_science_spw
   - total_science_bw_hz
   - min_channel_width_hz
   - median_channel_width_hz
   - max_channel_width_hz
   - min_velocity_resolution_kms
   - median_velocity_resolution_kms
   - max_velocity_resolution_kms
   - min_spw_bandwidth_hz
   - median_spw_bandwidth_hz
   - max_spw_bandwidth_hz
   - frequency_span_hz
   - mean_line_range_fraction_of_spw across SPWs
   - mean_buffered_line_range_fraction_of_spw across SPWs
   - mean_summed_line_fraction_of_spw across SPWs
   - mean_buffered_summed_line_fraction_of_spw across SPWs
   - number of distinct SPW bandwidth classes
   - number of distinct velocity-resolution classes
   - canonical setup signature
   - setup family
   - EB consistency flag

5. Summary/statistics outputs:
   - `summary_overall.txt` or `summary_overall.md`
   - `signature_summary_overall.csv`
   - `signature_summary_by_science_category.csv`, if metadata CSV is supplied
   - `signature_summary_by_band.csv`, if metadata CSV is supplied
   - `signature_summary_by_array_or_config.csv`, if metadata CSV is supplied
   - `parse_errors.csv`

Binning and signatures:
1. Bin each SPW by:
   - velocity-resolution bin: fine / medium / coarse
   - SPW-total-bandwidth bin: narrow / medium / wide
2. Use either data-driven tertiles or supplied physical bin edges. Make this configurable.
3. Assign each SPW a class such as:
   - fine-narrow
   - fine-medium
   - fine-wide
   - medium-narrow
   - medium-medium
   - medium-wide
   - coarse-narrow
   - coarse-medium
   - coarse-wide
4. For each MOUS setup, create an unordered canonical signature, e.g.:
   - `4x_coarse-wide`
   - `3x_coarse-wide + 1x_fine-narrow`
   - `2x_coarse-wide + 2x_fine-narrow`
5. Do not assume SPW order is meaningful.
6. Do not assume NSPW = 4. Measure the full NSPW distribution.

Optional broad setup-family labels:
1. Derive a broad `setup_family` from the SPW classes and NSPW:
   - continuum_like
   - mixed_line_continuum
   - line_heavy
   - many_spw_complex
   - full_pol, if polarization information supports this
   - other_or_unclassified
2. Keep this rule-based and transparent. Put the rules in the output summary.

Sanity checks and high-level statistics:
1. Report total:
   - number of MOUS directories found
   - number of listobs files found
   - number of cont.dat files found
   - number of EBs parsed
   - number of SPWs parsed
   - number/fraction of MOUS with usable listobs data
   - number/fraction of MOUS with usable cont.dat data
2. Report NSPW distribution:
   - count and fraction of MOUS with 1, 2, 3, 4, 5, 6, >6 science SPWs
3. Report EB consistency within MOUS:
   - fraction of MOUS where all EB spectral setups match
   - fraction where EBs differ in SPW count
   - fraction where EBs differ in nchan, bandwidth, or channel width
   - examples of inconsistent MOUS
4. Report distributions:
   - nchan
   - total_bw_hz
   - channel_width_hz
   - velocity_resolution_kms
   - line_range_fraction_of_spw
   - buffered_line_range_fraction_of_spw
   - summed_line_fraction_of_spw
   - buffered_summed_line_fraction_of_spw
5. Report top 5–10 setup signatures overall:
   - count
   - fraction
   - cumulative fraction
6. If metadata CSV is supplied, report top signatures by:
   - science category
   - ALMA band
   - array type / config proxy
7. Report parse failures and ambiguous cases clearly.
8. Include a few representative parsed examples in the human-readable summary for manual inspection.

Important assumptions:
1. Treat `listobs` as the authoritative source for native visibility spectral metadata.
2. Treat `cont.dat` as an annotation for continuum-vs-line coverage, not as the source of native SPW channelization.
3. Treat non-continuum ranges, i.e. the complement of cont.dat continuum ranges within an SPW, as line ranges.
4. Compute both:
   - span-based line coverage: min line frequency to max line frequency;
   - summed-disjoint line coverage: total width of all non-continuum intervals.
5. Also compute 10-channel buffered versions of both line-coverage metrics.
6. For a MOUS-level line-range result, average over targets and EBs where applicable.
7. If EBs within a MOUS have consistent spectral setups, the first EB can be treated as representative for setup classification.
8. Preserve the EB-level data anyway so consistency can be checked and future analyses can revisit the aggregation.
9. Output should be both machine-readable and human-readable.
10. Prefer robust parsing and clear warnings over silent assumptions.
```

