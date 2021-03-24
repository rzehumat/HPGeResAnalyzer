import pandas as pd
import glob
import numpy as np
import os
import re

from math import modf, isnan
from pathlib import Path
from sys import argv

def file_name_parse(file_name):
    name = file_name.split("/")[-1].split(".")[0]
    isotope, geometry = name.split("_g")
    A, element = re.split(r'(\d+)(\w+)', isotope)#[-3:-1]
    return A, element, geometry

def add_epsilon_dir(parsed_dir):
    if not os.path.isdir(parsed_dir):
        raise Exception(f"Directory '{parsed_dir}' not found. Consider creating it and moving parsed files there.")
    
    CALIBRATION_PATH = "aux_data/epsilons.csv"
    eps_df = pd.read_csv(CALIBRATION_PATH, index_col=0)
    Path("./with_eps").mkdir(parents=True, exist_ok=True)
    for parsed_file in glob.iglob(f"{parsed_dir}/*.csv"):
        geom = (parsed_file.split('.')[-2]).split('_g')[-1]
        print(geom)
        print(f"Geom is {geom} mm.")
        add_epsilon_file(parsed_file, geom, eps_df)

def drop_unit(long_str):
    return float(long_str.split(' ')[0])

# def add_epsilon_file(report_path, geom, eps_df):
def add_epsilon_file(df, geom, eps_df):
    # df = pd.read_csv(report_path, index_col=0)
    # file_name = report_path.split('/')[1]
    try:
        df["eps"] = df["E_tab"].apply(lambda x: add_epsilon_val(x, geom, eps_df))
    except:
        print("Seems like this isotope has no gamma lines. WTF have you measured.")
        print("I'll at least add detection efficiencies...")
        df["E_tab"] = df["Energy"]
        df["eps"] = df["E_tab"].apply(lambda x: add_epsilon_val(x, geom, eps_df))

    # change the order of columns
    # UGLY
    if "Ig" not in df.columns:
        df["Ig"] = float("NaN")
    
    # df["Live Time"] = df["Live Time"].apply(drop_unit)
    # df["Real Time"] = df["Real Time"].apply(drop_unit)
    # df["Dead Time"] = 0.001*df["Dead Time"].apply(drop_unit)
    # df["Identification Energy Tolerance"] = df["Identification Energy Tolerance"].apply(drop_unit)


    # old_column_order = df.columns.tolist()
    df["Geometry"] = geom
    # iso = (report_path.split('/')[-1]).split('_g')[0]
    # df["Isotope"] = iso
    # first_cols = ["Pk", "Isotope", "Energy", "FWHM", "E_tab", "Area", "%err", "Ig", "eps", "Geometry", "Real Time", "Live Time"]
    
    # new_cols = first_cols + list(set(old_column_order) - set(first_cols))
    # df = df[new_cols]
    # that is wrong -- it's always 30 in the RPT files...
    # to_drop = [
    #     "Sample Geometry", 
    #     "Use Fixed FWHM", 
    #     "Peak Analysis Report                    26.11.2020  5", 
    #     "Max Iterations", 
    #     "Peak Search Sensitivity", 
    #     "Peak Analysis From Channel", 
    #     "Peak Fit Engine Name",
    #     "Sample Size",
    #     "Sample Type", 
    #     "Sample Identification", 
    #     "IT", 
    #     "Fit"
    #     ]
    # df = df.drop(columns = to_drop)

    print("Epsilons added")

    return df

    # df.to_csv(f"out/{file_name}")

def add_epsilon_val(num, geom, eps_df):
    if isnan(num):
        return float("NaN")
    else:
        dec_e, int_e = modf(num)
        int_e = int(int_e)
        
        curr_eps = float(eps_df.loc[int_e, f"eps_{geom}_mm"])
        eps = curr_eps + dec_e*(float(eps_df.loc[int_e + 1, f"eps_{geom}_mm"]) - curr_eps)
    return eps

if __name__ == "__main__":
    file_dir = argv[1]
    add_epsilon_dir(file_dir)