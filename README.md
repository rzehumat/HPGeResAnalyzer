# HPGeResAnalyzer

![Lines of code](https://img.shields.io/tokei/lines/github/rzehumat/HPGeResAnalyzer)

A minimal script to parse Genie2k RPT files and process measured spectrum

## TLDR;
- insert `RPT` files to `raw_reports` dir
- run `python main.py` and add follow the instructions
- find results in `out` dir or report an issue

## Scripts
- `rptParser.py` parses data from `.RPT` file into a `csv`
- `getIg.py` adds tablular energy values and intensities of gamma-lines to individual peaks (if isotope is known -- unknown isotopes are TODO)
  - if there is no `csv` with gammas from given isotope, it is downloaded from _yellow pages_
