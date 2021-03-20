# HPGeResAnalyzer
A minimal script to parse Genie2k RPT files and process measured spectrum.

## Scripts
- `rptParser.py` parses data from `.RPT` file into a `csv`
- `getIg.py` adds tablular energy values and intensities of gamma-lines to individual peaks (if isotope is known -- unknown isotopes are TODO)
  - if there is no `csv` with gammas from given isotope, it is downloaded from _yellow pages_
