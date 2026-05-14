## EB-Level Volume Fractions

| Bucket | Mean | Median | Q25 | Q75 | IQR |
| --- | --- | --- | --- | --- | --- |
| science | 0.6847 | 0.6778 | 0.5236 | 0.9295 | 0.4059 |
| phase | 0.0728 | 0.0711 | 0.0191 | 0.1009 | 0.0818 |
| bandpass | 0.1155 | 0.0880 | 0.0302 | 0.1515 | 0.1214 |
| check_source | 0.0104 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| pointing | 0.0176 | 0.0073 | 0.0014 | 0.0207 | 0.0193 |
| polarization | 0.0012 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| diffgain_reference | 0.0002 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| diffgain_on_source | 0.0005 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| atmosphere | 0.0970 | 0.0451 | 0.0108 | 0.1462 | 0.1354 |
| other | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |


# Spectral Setup Signature Summary

## Inputs

- Unique EB-SPW rows used for binning: 28253
- Total MOUS assigned a setup: 2336
- Total EBs represented by those MOUS: 5559
- Total EB-SPWs represented by those EBs: 28253

## Cut Definitions

- Bandwidth cuts (GHz, fixed): low=0.090000, high=1.200000
- Resolution cuts (km/s): q33=0.593903, q66=7.736921, applied at 0.593904 between 0.593900 and 0.593908, and 7.741309 between 7.736921 and 7.745696
- Token convention: first letter is bandwidth (`N`, `M`, `W`), second letter is resolution (`F`, `M`, `C`).
- EB data rate uses fiducial Nant = 43 for 12m and 10 for 7m, with 4 bytes per complex correlation product and the science-scan median integration time per SPW.

## Category Counts

- Bandwidth bins over EB-SPW rows: {'medium': 5622, 'narrow': 5927, 'wide': 16704}
- Resolution bins over EB-SPW rows: {'coarse': 9416, 'fine': 9418, 'mid': 9419}

## Top 34 Setup Signatures (covering 85% of EBs target)

- [WC,WC,WC,WC]: 814 MOUS (34.85%), 1494 EBs (26.88%), 5976 EB-SPWs (21.15%), mean EB data rate = 0.0009 GB/s, line fraction mean/median/stddev = 22.90% / 18.75% / 11.93%
- [WM,WM,WM,WM]: 395 MOUS (16.91%), 869 EBs (15.63%), 3476 EB-SPWs (12.30%), mean EB data rate = 0.0046 GB/s, line fraction mean/median/stddev = 35.85% / 32.66% / 21.45%
- [MM,MM,MM,WC,WM]: 33 MOUS (1.41%), 209 EBs (3.76%), 1045 EB-SPWs (3.70%), mean EB data rate = 0.0023 GB/s, line fraction mean/median/stddev = 25.04% / 23.51% / 7.60%
- [WC,WC,WC,WM]: 69 MOUS (2.95%), 172 EBs (3.09%), 688 EB-SPWs (2.44%), mean EB data rate = 0.0024 GB/s, line fraction mean/median/stddev = 20.27% / 19.43% / 8.18%
- [MM,WC,WC,WC]: 32 MOUS (1.37%), 160 EBs (2.88%), 640 EB-SPWs (2.27%), mean EB data rate = 0.0012 GB/s, line fraction mean/median/stddev = 20.73% / 18.98% / 6.60%
- [MC,WC,WC,WC]: 52 MOUS (2.23%), 147 EBs (2.64%), 588 EB-SPWs (2.08%), mean EB data rate = 0.0017 GB/s, line fraction mean/median/stddev = 20.46% / 16.72% / 11.59%
- [WC,WM,WM,WM]: 34 MOUS (1.46%), 134 EBs (2.41%), 536 EB-SPWs (1.90%), mean EB data rate = 0.0027 GB/s, line fraction mean/median/stddev = 26.86% / 26.53% / 12.59%
- [MF,NF,NF,NF,NF,NF,NF,NF,NF,WM]: 13 MOUS (0.56%), 124 EBs (2.23%), 1240 EB-SPWs (4.39%), mean EB data rate = 0.0030 GB/s, line fraction mean/median/stddev = 41.65% / 28
.96% / 27.98%
- [MM,MM,WM,WM]: 28 MOUS (1.20%), 123 EBs (2.21%), 492 EB-SPWs (1.74%), mean EB data rate = 0.0004 GB/s, line fraction mean/median/stddev = 23.36% / 21.28% / 14.51%
- [MF,NF,NF,NF,NF,WM,WM]: 18 MOUS (0.77%), 115 EBs (2.07%), 805 EB-SPWs (2.85%), mean EB data rate = 0.0053 GB/s, line fraction mean/median/stddev = 51.45% / 55.55% / 13.
42%
- [WC,WC,WM,WM]: 63 MOUS (2.70%), 100 EBs (1.80%), 400 EB-SPWs (1.42%), mean EB data rate = 0.0042 GB/s, line fraction mean/median/stddev = 32.78% / 32.29% / 13.78%
- [NF,NF,NF,NF,NF,NF,NF,WM]: 26 MOUS (1.11%), 95 EBs (1.71%), 760 EB-SPWs (2.69%), mean EB data rate = 0.0073 GB/s, line fraction mean/median/stddev = 44.87% / 51.61% / 20.98%
- [MF,MF,WC,WC]: 26 MOUS (1.11%), 90 EBs (1.62%), 360 EB-SPWs (1.27%), mean EB data rate = 0.0013 GB/s, line fraction mean/median/stddev = 23.68% / 21.27% / 6.06%
- [MF,MF,NF,NF,NF,NF,NF,NF,WM]: 19 MOUS (0.81%), 76 EBs (1.37%), 684 EB-SPWs (2.42%), mean EB data rate = 0.0005 GB/s, line fraction mean/median/stddev = 19.63% / 17.09% / 11.96%
- [MF,MM,MM,WC]: 6 MOUS (0.26%), 73 EBs (1.31%), 292 EB-SPWs (1.03%), mean EB data rate = 0.0025 GB/s, line fraction mean/median/stddev = 16.80% / 16.40% / 2.47%
- [WM,WM,WM,WM,WM,WM,WM,WM]: 39 MOUS (1.67%), 71 EBs (1.28%), 568 EB-SPWs (2.01%), mean EB data rate = 0.0025 GB/s, line fraction mean/median/stddev = 32.78% / 33.20% / 14.28%
- [MC,MC,WC,WC]: 29 MOUS (1.24%), 63 EBs (1.13%), 252 EB-SPWs (0.89%), mean EB data rate = 0.0008 GB/s, line fraction mean/median/stddev = 27.55% / 26.45% / 6.81%
- [MF,WM,WM,WM]: 36 MOUS (1.54%), 60 EBs (1.08%), 240 EB-SPWs (0.85%), mean EB data rate = 0.0060 GB/s, line fraction mean/median/stddev = 16.80% / 11.87% / 18.43%
- [MM,NF,NF,NF]: 9 MOUS (0.39%), 46 EBs (0.83%), 184 EB-SPWs (0.65%), mean EB data rate = 0.0093 GB/s, line fraction mean/median/stddev = 30.05% / 32.00% / 13.68%
- [MF,MF,MF,MF,WM]: 17 MOUS (0.73%), 44 EBs (0.79%), 220 EB-SPWs (0.78%), mean EB data rate = 0.0005 GB/s, line fraction mean/median/stddev = 31.25% / 19.77% / 20.95%
- [NF,NF,NF,NF,NF,NF,NF,NF,WM]: 16 MOUS (0.68%), 42 EBs (0.76%), 378 EB-SPWs (1.34%), mean EB data rate = 0.0063 GB/s, line fraction mean/median/stddev = 27.55% / 24.45% / 12.58%
- [MM,MM,WC,WC,WM]: 29 MOUS (1.24%), 42 EBs (0.76%), 210 EB-SPWs (0.74%), mean EB data rate = 0.0049 GB/s, line fraction mean/median/stddev = 12.41% / 12.30% / 2.24%
- [MF,NF,NF,NF,NF,NF,NF,NF,NF,NF,NF,WM]: 26 MOUS (1.11%), 39 EBs (0.70%), 468 EB-SPWs (1.66%), mean EB data rate = 0.0092 GB/s, line fraction mean/median/stddev = 31.13% / 26.68% / 13.58%
- [MF,MF,MF,MF,MF,MF,MF,MF,MF,MF,MF,MF,MF,MF]: 19 MOUS (0.81%), 38 EBs (0.68%), 532 EB-SPWs (1.88%), mean EB data rate = 0.0009 GB/s, line fraction mean/median/stddev = 26.84% / 24.69% / 7.16%
- [NF,NF,NF,NF,NF,NF,NF,NF,NF,NF,NF,NF,WM]: 24 MOUS (1.03%), 38 EBs (0.68%), 494 EB-SPWs (1.75%), mean EB data rate = 0.0097 GB/s, line fraction mean/median/stddev = 29.49% / 28.83% / 11.17%
- [MF,MF,MF,WM]: 15 MOUS (0.64%), 37 EBs (0.67%), 148 EB-SPWs (0.52%), mean EB data rate = 0.0107 GB/s, line fraction mean/median/stddev = 33.74% / 35.03% / 16.54%
- [MF,MF,MF,MF,MF,MF,MF,WM]: 11 MOUS (0.47%), 34 EBs (0.61%), 272 EB-SPWs (0.96%), mean EB data rate = 0.0005 GB/s, line fraction mean/median/stddev = 22.36% / 16.10% / 12.87%
- [MM,WM,WM,WM]: 17 MOUS (0.73%), 32 EBs (0.58%), 128 EB-SPWs (0.45%), mean EB data rate = 0.0089 GB/s, line fraction mean/median/stddev = 37.72% / 35.11% / 16.21%
- [WM,WM,WM,WM,WM]: 3 MOUS (0.13%), 32 EBs (0.58%), 128 EB-SPWs (0.45%), mean EB data rate = 0.0010 GB/s, line fraction mean/median/stddev = nan% / nan% / nan%
- [MF,MF,MF,MF,MM,WM]: 15 MOUS (0.64%), 29 EBs (0.52%), 174 EB-SPWs (0.62%), mean EB data rate = 0.0005 GB/s, line fraction mean/median/stddev = 8.50% / 8.07% / 3.80%
- [MF,MF,WM,WM,WM]: 14 MOUS (0.60%), 28 EBs (0.50%), 140 EB-SPWs (0.50%), mean EB data rate = 0.0045 GB/s, line fraction mean/median/stddev = 19.77% / 13.47% / 17.73%
- [MF,MF,MF,MF]: 19 MOUS (0.81%), 28 EBs (0.50%), 112 EB-SPWs (0.40%), mean EB data rate = 0.0216 GB/s, line fraction mean/median/stddev = 50.65% / 53.01% / 26.68%
- [NF,NF,NF,WM]: 13 MOUS (0.56%), 28 EBs (0.50%), 112 EB-SPWs (0.40%), mean EB data rate = 0.0005 GB/s, line fraction mean/median/stddev = 11.03% / 12.77% / 3.50%
- [NF,NF,NF,NF,NF,NF,NF,NF,NF,NF,NF,NF,WF]: 14 MOUS (0.60%), 27 EBs (0.49%), 351 EB-SPWs (1.24%), mean EB data rate = 0.0238 GB/s, line fraction mean/median/stddev = 31.92% / 22.12% / 21.44%

# ryanCy11 WSU Projection Summary

## Methods

- `memo_uniform_binned`: memo procedure using finest current SPW per MOUS and the 5 stepped2 bins.
- `memo_uniform_exact`: same, but keep the current finest requested resolution instead of snapping to the bin floor.
- `distributed_binned`: preserve the current within-MOUS SPW resolution distribution, but map each SPW into the 5 stepped2 bins.
- `distributed_exact`: preserve the current within-MOUS SPW resolution distribution using the current requested resolution directly.

## Overall Results

### M1

- memo_uniform_binned: sample-total factor 18.431x, mean EB factor 9.825x, median EB factor 4.800x
- memo_uniform_exact: sample-total factor 16.754x, mean EB factor 7.345x, median EB factor 2.178x
- distributed_binned: sample-total factor 10.705x, mean EB factor 6.508x, median EB factor 3.964x
- distributed_exact: sample-total factor 8.993x, mean EB factor 4.316x, median EB factor 1.582x

- Change from keeping exact finest resolution instead of bin-flooring: 0.909x relative to memo_uniform_binned.
- Change from preserving the SPW resolution distribution at binned resolution: 0.581x relative to memo_uniform_binned.
- Additional change from using exact per-SPW resolutions on top of preserving the distribution: 0.840x relative to distributed_binned.

### M4

- memo_uniform_binned: sample-total factor 33.780x, mean EB factor 15.795x, median EB factor 7.375x
- memo_uniform_exact: sample-total factor 30.469x, mean EB factor 10.523x, median EB factor 2.979x
- distributed_binned: sample-total factor 20.588x, mean EB factor 11.156x, median EB factor 5.276x
- distributed_exact: sample-total factor 16.679x, mean EB factor 6.389x, median EB factor 2.347x

- Change from keeping exact finest resolution instead of bin-flooring: 0.902x relative to memo_uniform_binned.
- Change from preserving the SPW resolution distribution at binned resolution: 0.609x relative to memo_uniform_binned.
- Additional change from using exact per-SPW resolutions on top of preserving the distribution: 0.810x relative to distributed_binned.

### M5

- memo_uniform_binned: sample-total factor 89.974x, mean EB factor 38.573x, median EB factor 14.375x
- memo_uniform_exact: sample-total factor 82.014x, mean EB factor 27.153x, median EB factor 8.517x
- distributed_binned: sample-total factor 53.163x, mean EB factor 26.262x, median EB factor 11.049x
- distributed_exact: sample-total factor 43.974x, mean EB factor 16.114x, median EB factor 5.078x

- Change from keeping exact finest resolution instead of bin-flooring: 0.912x relative to memo_uniform_binned.
- Change from preserving the SPW resolution distribution at binned resolution: 0.591x relative to memo_uniform_binned.
- Additional change from using exact per-SPW resolutions on top of preserving the distribution: 0.827x relative to distributed_binned.

