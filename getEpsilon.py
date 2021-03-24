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

def add_epsilon_file(df, geom, eps_df):
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
    
    df["Geometry"] = geom
    print("Epsilons added")

    return df


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