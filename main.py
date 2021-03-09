import rptParser
import getIg
import getEpsilon
#import searchRadiation

print("Available modes:")
print("1 ... 'Nuclide search': We know the isotope, peak detection efficiency and geometry")
print("Parses RPT files, downloads and adds I_gamma to peaks, adds epsilons to peaks and saves as csv.")
print("2 ... 'Radiation search': We do NOT know the isotopes, peak detection efficiency and geometry")
print("Parses RPT files, searches the radiation, adds possible isotopes to peaks, adds I_gamma to peaks, adds epsilons to peaks and saves as csv.")

print("Select mode:")
mode = input("[1] ")

if mode == "1":
    print("Expected filenames are like 1H_g80.RPT")
    # rpt_dir = input("Relative path to directory with RPT files: ")
    rptParser.parse_RPT("raw_reports")
    getIg.append_Igamma_dir("parsed_reports/")
    getEpsilon.add_epsilon_dir("with_Ig/")
elif mode == 2:
    rptParser.parse_RPT("naa_reports")
    searchRadiation.search_dir("parsed_reports")