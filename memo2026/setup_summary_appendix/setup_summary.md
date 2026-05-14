# Spectral Setup Signature Summary

## Inputs

- Unique EB-SPW rows used for binning: 28253
- CONT EB-SPW rows excluded before deriving non-CONT cuts: 6107
- Non-CONT EB-SPW rows used to derive cuts: 22146
- Total MOUS assigned a setup: 2336
- Total EBs represented by those MOUS: 5559
- Total EB-SPWs represented by those EBs: 28253

## Cut Definitions

- Bandwidth cuts (GHz, fixed): low=0.09, high=1.2
- CONT definition: nchan <= 128 and bandwidth >= 1.88 GHz.
- Resolution cuts (km/s): q33=0.317, q66=1.81, applied at 0.317 between 0.317 and 0.317, and 1.81 between 1.81 and 1.81
- Token convention: `CONT` marks TDM-like continuum SPWs; non-CONT tokens use bandwidth (`N`, `M`, `W`) plus resolution (`F`, `M`, `C`).
- EB data rate uses fiducial Nant = 43 for 12m and 10 for 7m, with 4 bytes per complex correlation product and the science-scan median integration time per SPW.

## Category Counts

- CONT count over EB-SPW rows: 6107
- Bandwidth bins over non-CONT EB-SPW rows: {'medium': 5622, 'narrow': 5927, 'wide': 10597}
- Resolution bins over non-CONT EB-SPW rows: {'coarse': 7382, 'fine': 7382, 'mid': 7382}

## Mixed Resolution At EB Level

- COARSE: 1458 EBs (26.2%), single-class
- CONT: 851 EBs (15.3%), single-class
- MEDIUM: 591 EBs (10.6%), single-class
- FINE: 73 EBs (1.31%), single-class
- FINE + MEDIUM: 1083 EBs (19.5%), mixed
- CONT + COARSE: 550 EBs (9.89%), mixed
- CONT + MEDIUM: 546 EBs (9.82%), mixed
- MEDIUM + COARSE: 207 EBs (3.72%), mixed
- CONT + FINE: 88 EBs (1.58%), mixed
- FINE + COARSE: 85 EBs (1.53%), mixed
- CONT + MEDIUM + COARSE: 14 EBs (0.252%), mixed
- CONT + FINE + COARSE: 10 EBs (0.18%), mixed
- CONT + FINE + MEDIUM: 2 EBs (0.036%), mixed
- FINE + MEDIUM + COARSE: 1 EBs (0.018%), mixed

## Top 40 Setup Signatures (covering 85% of EBs target)

- [WC,WC,WC,WC]: 394 MOUS (16.9%), 1047 EBs (18.8%), 4188 EB-SPWs (14.8%), mean EB data rate = 0.0023 GB/s, fully empty fraction = 0%, line fraction mean/median/stddev = 33.6% / 32.8% / 18.2%
- [CONT,CONT,CONT,CONT]: 511 MOUS (21.9%), 814 EBs (14.6%), 3256 EB-SPWs (11.5%), mean EB data rate = 0.00098 GB/s, fully empty fraction = 0%, line fraction mean/median/stddev = 21.4% / 17.8% / 10.1%
- [WM,WM,WM,WM]: 221 MOUS (9.46%), 421 EBs (7.57%), 1684 EB-SPWs (5.96%), mean EB data rate = 0.00688 GB/s, fully empty fraction = 0%, line fraction mean/median/stddev = 31.5% / 20.9% / 23.1%
- [CONT,CONT,CONT,WC]: 147 MOUS (6.29%), 228 EBs (4.1%), 912 EB-SPWs (3.23%), mean EB data rate = 0.00159 GB/s, fully empty fraction = 0%, line fraction mean/median/stddev = 20.1% / 19.7% / 5.21%
- [CONT,MM,MM,MM,WM]: 33 MOUS (1.41%), 209 EBs (3.76%), 1045 EB-SPWs (3.7%), mean EB data rate = 0.0023 GB/s, fully empty fraction = 0%, line fraction mean/median/stddev = 25% / 23.5% / 7.6%
- [MC,WC,WC,WC]: 48 MOUS (2.05%), 164 EBs (2.95%), 656 EB-SPWs (2.32%), mean EB data rate = 0.00253 GB/s, fully empty fraction = 0%, line fraction mean/median/stddev = 20.7% / 16% / 11.5%
- [CONT,CONT,CONT,MC]: 40 MOUS (1.71%), 154 EBs (2.77%), 616 EB-SPWs (2.18%), mean EB data rate = 0.000256 GB/s, fully empty fraction = 0%, line fraction mean/median/stddev = 21% / 18.7% / 8.16%
- [MF,NF,NF,NF,NF,NF,NF,NF,NF,WM]: 13 MOUS (0.557%), 124 EBs (2.23%), 1240 EB-SPWs (4.39%), mean EB data rate = 0.00299 GB/s, fully empty fraction = 0%, line fraction mean/median/stddev = 41.6% / 29% / 28%
- [MF,NF,NF,NF,NF,WM,WM]: 18 MOUS (0.771%), 115 EBs (2.07%), 805 EB-SPWs (2.85%), mean EB data rate = 0.00531 GB/s, fully empty fraction = 0%, line fraction mean/median/stddev = 51.4% / 55.5% / 13.4%
- [MC,MC,WC,WC]: 15 MOUS (0.642%), 111 EBs (2%), 444 EB-SPWs (1.57%), mean EB data rate = 0.000676 GB/s, fully empty fraction = 0%, line fraction mean/median/stddev = 23.7% / 22.4% / 15.2%
- [NF,NF,NF,NF,NF,NF,NF,WM]: 26 MOUS (1.11%), 95 EBs (1.71%), 760 EB-SPWs (2.69%), mean EB data rate = 0.00727 GB/s, fully empty fraction = 0%, line fraction mean/median/stddev = 44.9% / 51.6% / 21%
- [CONT,CONT,CONT,WM]: 21 MOUS (0.899%), 86 EBs (1.55%), 344 EB-SPWs (1.22%), mean EB data rate = 0.000544 GB/s, fully empty fraction = 0%, line fraction mean/median/stddev = 17% / 14.6% / 6.61%
- [MF,MF,NF,NF,NF,NF,NF,NF,WM]: 19 MOUS (0.813%), 76 EBs (1.37%), 684 EB-SPWs (2.42%), mean EB data rate = 0.000457 GB/s, fully empty fraction = 0%, line fraction mean/median/stddev = 19.6% / 17.1% / 12%
- [CONT,CONT,MM,MM]: 21 MOUS (0.899%), 73 EBs (1.31%), 292 EB-SPWs (1.03%), mean EB data rate = 7.47e-05 GB/s, fully empty fraction = 0%, line fraction mean/median/stddev = 24% / 21.4% / 5.17%
- [CONT,MM,MM,MM]: 6 MOUS (0.257%), 73 EBs (1.31%), 292 EB-SPWs (1.03%), mean EB data rate = 0.00246 GB/s, fully empty fraction = 0%, line fraction mean/median/stddev = 16.8% / 16.4% / 2.47%
- [CONT,CONT,WC,WC]: 34 MOUS (1.46%), 66 EBs (1.19%), 264 EB-SPWs (0.934%), mean EB data rate = 0.00104 GB/s, fully empty fraction = 0%, line fraction mean/median/stddev = 25.7% / 25.1% / 8.12%
- [NF,NF,NF,NF,NF,NF,NF,NF,NF,NF,NF,NF,WM]: 38 MOUS (1.63%), 65 EBs (1.17%), 845 EB-SPWs (2.99%), mean EB data rate = 0.0155 GB/s, fully empty fraction = 0%, line fraction mean/median/stddev = 30.4% / 26.3% / 15.8%
- [CONT,CONT,MC,MC]: 29 MOUS (1.24%), 65 EBs (1.17%), 260 EB-SPWs (0.92%), mean EB data rate = 0.000875 GB/s, fully empty fraction = 0%, line fraction mean/median/stddev = 28.7% / 27.1% / 6.57%
- [MF,WM,WM,WM]: 37 MOUS (1.58%), 60 EBs (1.08%), 240 EB-SPWs (0.849%), mean EB data rate = 0.00774 GB/s, fully empty fraction = 0%, line fraction mean/median/stddev = 24.2% / 12.2% / 26.7%
- [NF,NF,NF,NF,NF,NF,NF,NF,WM]: 22 MOUS (0.942%), 59 EBs (1.06%), 531 EB-SPWs (1.88%), mean EB data rate = 0.0113 GB/s, fully empty fraction = 0%, line fraction mean/median/stddev = 26.6% / 20.3% / 13.8%
- [WC,WC,WC,WM]: 22 MOUS (0.942%), 52 EBs (0.935%), 208 EB-SPWs (0.736%), mean EB data rate = 0.00143 GB/s, fully empty fraction = 0%, line fraction mean/median/stddev = 37.7% / 42.4% / 14.4%
- [MF,MF,MF,MF,WM]: 17 MOUS (0.728%), 44 EBs (0.792%), 220 EB-SPWs (0.779%), mean EB data rate = 0.000457 GB/s, fully empty fraction = 0%, line fraction mean/median/stddev = 31.2% / 19.8% / 20.9%
- [CONT,CONT,MM,MM,WM]: 29 MOUS (1.24%), 42 EBs (0.756%), 210 EB-SPWs (0.743%), mean EB data rate = 0.00489 GB/s, fully empty fraction = 0%, line fraction mean/median/stddev = 12.4% / 12.3% / 2.24%
- [MF,NF,NF,NF,NF,NF,NF,NF,NF,NF,NF,WM]: 26 MOUS (1.11%), 39 EBs (0.702%), 468 EB-SPWs (1.66%), mean EB data rate = 0.00917 GB/s, fully empty fraction = 0%, line fraction mean/median/stddev = 31.1% / 26.7% / 13.6%
- [CONT,CONT,CONT,CONT,CONT,CONT,CONT,CONT]: 21 MOUS (0.899%), 36 EBs (0.648%), 276 EB-SPWs (0.977%), mean EB data rate = 0.00297 GB/s, fully empty fraction = 0%, line fraction mean/median/stddev = 19.3% / 18.6% / 4.11%
- [MM,NF,NF,NF]: 6 MOUS (0.257%), 36 EBs (0.648%), 144 EB-SPWs (0.51%), mean EB data rate = 0.0101 GB/s, fully empty fraction = 0%, line fraction mean/median/stddev = 31.2% / 36.6% / 13.4%
- [MM,MM,MM,MM,MM,MM,MM,WM]: 11 MOUS (0.471%), 34 EBs (0.612%), 272 EB-SPWs (0.963%), mean EB data rate = 0.000457 GB/s, fully empty fraction = 0%, line fraction mean/median/stddev = 22.4% / 16.1% / 12.9%
- [WC,WC,WC,WC,WC]: 3 MOUS (0.128%), 34 EBs (0.612%), 136 EB-SPWs (0.481%), mean EB data rate = 0.000543 GB/s, fully empty fraction = 0%, line fraction mean/median/stddev = nan% / nan% / nan%
- [MF,MF,MF,MF,MF,MM,MM,MM,MM,MM,MM,MM,MM,MM]: 16 MOUS (0.685%), 32 EBs (0.576%), 448 EB-SPWs (1.59%), mean EB data rate = 0.000914 GB/s, fully empty fraction = 0%, line fraction mean/median/stddev = 27.5% / 24.9% / 7.29%
- [MM,MM,WC,WC]: 19 MOUS (0.813%), 31 EBs (0.558%), 124 EB-SPWs (0.439%), mean EB data rate = 0.000457 GB/s, fully empty fraction = 0%, line fraction mean/median/stddev = 27.6% / 26% / 15.5%
- [MF,MF,WM,WM,WM]: 16 MOUS (0.685%), 30 EBs (0.54%), 150 EB-SPWs (0.531%), mean EB data rate = 0.00461 GB/s, fully empty fraction = 0%, line fraction mean/median/stddev = 13% / 12.4% / 4.22%
- [MM,MM,MM,MM,MM,WC]: 15 MOUS (0.642%), 29 EBs (0.522%), 174 EB-SPWs (0.616%), mean EB data rate = 0.000457 GB/s, fully empty fraction = 0%, line fraction mean/median/stddev = 8.5% / 8.07% / 3.8%
- [NF,NF,NF,WC]: 13 MOUS (0.557%), 28 EBs (0.504%), 112 EB-SPWs (0.396%), mean EB data rate = 0.000457 GB/s, fully empty fraction = 0%, line fraction mean/median/stddev = 11% / 12.8% / 3.5%
- [CONT,CONT,WM,WM]: 15 MOUS (0.642%), 27 EBs (0.486%), 108 EB-SPWs (0.382%), mean EB data rate = 0.00934 GB/s, fully empty fraction = 0%, line fraction mean/median/stddev = 40.6% / 36.6% / 11%
- [MF,MF,MF,MF]: 17 MOUS (0.728%), 26 EBs (0.468%), 104 EB-SPWs (0.368%), mean EB data rate = 0.0226 GB/s, fully empty fraction = 0%, line fraction mean/median/stddev = 47.7% / 44.9% / 26.7%
- [WC,WC,WC,WC,WC,WC,WC,WC]: 19 MOUS (0.813%), 24 EBs (0.432%), 192 EB-SPWs (0.68%), mean EB data rate = 0.00226 GB/s, fully empty fraction = 0%, line fraction mean/median/stddev = 42.5% / 38.5% / 8.8%
- [MM,NF,NF,NF,NF,NF]: 4 MOUS (0.171%), 24 EBs (0.432%), 144 EB-SPWs (0.51%), mean EB data rate = 0.000457 GB/s, fully empty fraction = 0%, line fraction mean/median/stddev = 23% / 22.7% / 3.09%
- [MM,MM,MM,WM,WM]: 13 MOUS (0.557%), 24 EBs (0.432%), 120 EB-SPWs (0.425%), mean EB data rate = 0.0198 GB/s, fully empty fraction = 0%, line fraction mean/median/stddev = 20.7% / 18.5% / 7.39%
- [MM,WM,WM,WM]: 13 MOUS (0.557%), 24 EBs (0.432%), 96 EB-SPWs (0.34%), mean EB data rate = 0.00882 GB/s, fully empty fraction = 0%, line fraction mean/median/stddev = 34.7% / 35.1% / 9.59%
- [MC,MC,MC,WC]: 11 MOUS (0.471%), 24 EBs (0.432%), 96 EB-SPWs (0.34%), mean EB data rate = 0.00171 GB/s, fully empty fraction = 0%, line fraction mean/median/stddev = 30.7% / 25.3% / 10.2%

## Coverage Of Top 40 Setups

- Top 40 setups cover 1999 / 2336 MOUS (85.6%).
- Top 40 setups cover 4745 / 5559 EBs (85.4%).
- Top 40 setups cover 23160 / 28253 EB-SPWs (82%).
