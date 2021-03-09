import pandas as pd
import glob
import numpy as np
import os
from math import modf, isnan
from pathlib import Path
from sys import argv

def add_epsilon_dir(parsed_dir, geom):
    if not os.path.isdir(parsed_dir):
        raise Exception(f"Directory '{parsed_dir}' not found. Consider creating it and moving parsed files there.")
    
    CALIBRATION_PATH = "aux_data/epsilons.csv"
    eps_df = pd.read_csv(CALIBRATION_PATH, index_col=0)
    Path("./with_eps").mkdir(parents=True, exist_ok=True)
    for parsed_file in glob.iglob(f"{parsed_dir}/*.csv"):
        add_epsilon_file(parsed_file, geom, eps_df)

def add_epsilon_file(report_path, geom, eps_df):
    df = pd.read_csv(report_path, index_col=0)
    file_name = report_path.split('/')[1]
    df["eps"] = df["E_tab"].apply(lambda x: add_epsilon_val(x, geom, eps_df))

    # change the order of columns 
    old_column_order = df.columns.tolist()
    df["Geometry"] = [geom for i in range(df.shape[0])] 
    first_cols = ["Pk", "Energy", "FWHM", "E_tab", "Area", "%err", "Ig", "eps", "Geometry"]
    new_cols = first_cols + list(set(old_column_order) - set(first_cols))
    df = df[new_cols]
    # that is wrong -- it's always 30 in the RPT files...
    df = df.drop(columns = ["Sample Geometry"])

    df.to_csv(f"with_eps/{file_name}")

def add_epsilon_val(num, geom, eps_df):
    if isnan(num):
        return float("NaN")
    else:
        dec_e, int_e = modf(num)
        int_e = int(int_e)
        
        curr_eps = float(eps_df.loc[int_e, f"eps_{geom}_mm"])
        eps = curr_eps + dec_e*(float(eps_df.loc[int_e + 1, f"eps_{geom}_mm"]) - curr_eps)
    return eps

file_dir = argv[1]
geom = argv[2]
add_epsilon_dir(file_dir, geom)