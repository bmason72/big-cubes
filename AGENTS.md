
This directory contains code, text, and latex files related to
projections of ALMA WSU data characteristics and-- using those
projections-- compute and infrastructure cost estimates.

Most directories contain an index.md file, briefly describing key
files in the directory, as well as secondary files. You should review
and comprehend the key files in relation to our tasks. The secondary
files can be looked at if additional context is needed. Other files
can also be looked at, but only if there is reason to consider that
they may contain additional, relevant information. Note:
  * You are free to add information to these index files. 
  * if there is a directory without an index file, you are free to make it.
  * note that for directories containing latex source, "main.tex" is usually
    the best place to start.

Our tasks:
1- highest level, I want a script that will correctly, defensibly, reproducibly and explainably
     calculate the data rates for M1, M4, and M5 -- as well as possible to-be-defined in the future
     stages; produce latex tables summarizing them; and produce plots illustrating their distribution.
     The tables and figures should also show the corresponding, current (BLC/ACA) quantities.
     The tables and figures in the WSU SDD (data properties appendix) are about the right level
     of detail and would be a good thing to aim to reproduce.
     We should also be able to reproduce the current BLC/ACA quantities.
     
2- I would like to know if there are any significant errors or
    inconsistencies in the projections or documents that we have not
    yet identified. The important documents are 1) the SDD data
    properties projections; 2) the TSI cost projection report; 3) the
    two WSU data properties memos.

3- I would like, in the near future, to make the projections more
   sophisticated.  Notable elements to include would be: "smart
   mitigation"; smarter channel-resolution assumptions (instead of
   assuming, say, all spectral windows are at the very highest
   resolution); smarter integration time choices (3sec for 12m long
   baseline; 6sec for not-long baseline probably C6 and shorter); using 3 instead
   of 5 pixels per synthesized beam.

4- I have considered that it could be helpful to have a GUI interface
   where we can explore the projections characteristics, potentially as
   a function of assumptions.  The python notebooks were essential that,
   but i find them fidgety and fragile. If there were a robust way to do this
   that does not add a lot of validation burden, it would be of value.
   It is the least important task.



