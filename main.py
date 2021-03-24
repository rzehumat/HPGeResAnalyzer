import rptParser
import getIg
import getEpsilon
import glob
import os

import pandas as pd

from pathlib import Path
#import searchRadiation

OUTPUT_DIR = "out"

print("Available modes (default 0):")
print("0 ... 'Process dir'")
print("1 ... 'Nuclide search': We know the isotope, peak detection efficiency and geometry")
print("Parses RPT files, downloads and adds I_gamma to peaks, adds epsilons to peaks and saves as csv.")
print("2 ... 'Radiation search': We do NOT know the isotopes, peak detection efficiency and geometry")
print("Parses RPT files, searches the radiation, adds possible isotopes to peaks, adds I_gamma to peaks, adds epsilons to peaks and saves as csv.")

print("Select mode:")
mode = input("[0]/[1]/[2] ") or "0"

if mode == "0":
    raw_dir = input("Enter relative path to directory (press enter to keep default 'raw_reports': ") or "raw_reports"
    if not os.path.isdir(raw_dir):
        raise Exception(f"Folder {raw_dir} not found. Consider creating it and moving RPT files there.")
    else:
        Path(OUTPUT_DIR).mkdir(parents=True, exist_ok=True)
        ig_all_df = pd.read_parquet("aux_data/ig_all.pq")
        CALIBRATION_PATH = "aux_data/epsilons.csv"
        eps_df = pd.read_csv(CALIBRATION_PATH, index_col=0)
        for raw_file in glob.iglob(f"{raw_dir}/*.RPT"):
            try:
                A, element, geometry = getEpsilon.file_name_parse(raw_file)
            except:
                print(f"Info not detected from file name '{raw_file.split('/')[-1]}'.")
                print("If known, add it manually. If not, leave blank.")
                A = input("A = ")
                element = input("Element (leave blank if unknown)= ")
                geometry = input("Geometry (3/30/80/120/250 mm): ")
            
            raw_df = rptParser.parse_one_RPT(raw_file)
            df_ig = getIg.append_Igamma(raw_df, A, element, ig_all_df)
            df_ig_eps = getEpsilon.add_epsilon_file(df_ig, geometry, eps_df)

            file_name = raw_file.split("/")[-1].split(".")[-2]
            df_ig_eps.to_csv(f"{OUTPUT_DIR}/{file_name}.csv", index=False)


# if mode == "1":
elif mode == "1":
    print("Expected filenames are like 1H_g80.RPT")
    # rpt_dir = input("Relative path to directory with RPT files: ")
    rptParser.parse_RPT("raw_reports")
    getIg.append_Igamma_dir("parsed_reports/")
    getEpsilon.add_epsilon_dir("with_Ig/")
elif mode == "2":
    rptParser.parse_RPT("naa_raw")
    #searchRadiation.search_dir("parsed_reports")