
Key files:
*sum_m1_prod.py - python script to calculate WSU Milestone 1 product sizes/characteristics

*Example_notebook_MS1_MS5.ipynb- python notebook used (I think) to calculate the numbers and plots
   regarding data characteristics in the WSU System Design Description.

*data - directory containing the database of ALMA project metadata which, extrapolated
   to WSU, are used to estimate the characteristics of the future datasets.

*wsu-sdd-excerpts - directory containing the SDD text, including the aforementioned data 
   characteristics that are used in the design and costing.

*Estimated_WSU_Data_Properties - directory containing latex source for the memo
   describing the methodology used to project WSU data characteristics.
   This is the original method memo, but the terminology and assumptions for the WSU
   different, future stages changed after it was written (but before the SDD).

*WSU_Estimated_Data_Properties_Update - directory containing latex source for the memo
   that provides updated data characteristics for the WSU project milestones used in the
   SDD (M1, M4, M5).
   Note that the estimated data properites (and preceding work, i think) *averaged* the 7m and 12m
    data rates, which is incorrect: they should be added.

*WSU\ SDD\ data\ volume\ -\ Stage\ Lookup.csv - table summarizing the key assumptions that differ
   between the SDD WSU scenarios (Milestone 1,4,5 = M1, M4, M5) and earlier terminology
   (IWS, MWS, FWS)

*README - original file documenting the data properties database and code

*24_10_03_DPIICRE_report.txt - Data processing cost and resource estimate TSI report.
   this report used the pre-SDD WSU milestone definitions, but is of interest since we would
   like in the near future to repeat this exercise - probably using a similar methodology
   to what was in the report - but using the new milestone definitions.
   (Cost was not in scope for the WSU PDR that the SDD was part of in 2025; near end of
   2026 we expect a project cost review)

Secondary files:

*various, earlier version of sum_m1_prod.py
  -sum_m1_prodByProj.py
  -sum_m1_prodGood.py
  -sum_m1_prodBetter.py

