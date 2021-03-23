import pandas as pd
import json
import numpy as np
import os
import glob
import re
import urllib

from pathlib import Path
from sys import argv
from getEpsilon import file_name_parse


def append_Igamma_dir(parsed_dir):
    if not os.path.isdir(parsed_dir):
        raise Exception(f"Directory '{parsed_dir}' not found. Consider creating it and moving parsed files there.")

    ig_df = pd.read_parquet("aux_data/ig_all.pq")
    for parsed_file in glob.iglob(f"{parsed_dir}/*.csv"):
        A, element, _ = file_name_parse(parsed_file)
        df = pd.read_csv(parsed_file, index_col=0)
        append_Igamma(df, A, element, ig_df)

def permute_columns(df, first_cols):
    old_column_order = df.columns.tolist()
    new_cols = first_cols + list(set(old_column_order) - set(first_cols))
    df = df[new_cols]
    return df

def append_Igamma(parsed_df, A, element, ig_all_df):
    # print(file)
    # A, element = re.split(r'(\d+)(\w+)', file)[-3:-1]
    # element = element.split("_")[0]
    # print(A)
    # print(element)

    if len(A) > 0:
        ig_df = ig_all_df.loc[f"{A}{element}"]
    else:
        ig_df = ig_all_df
    
    # parsed_df = pd.read_csv(file, index_col=0)
    joined_df = add_Ig(parsed_df, ig_df)
    
    Path("./with_Ig").mkdir(parents=True, exist_ok=True)
    joined_df = permute_columns(joined_df, ["Energy", "E_tab", "Ig", "Area", "sigm_E", "sigm_Ig", "FWHM", "%err", "Live Time", "Real Time", "Dead Time (rel)"])
    joined_df = joined_df.sort_values(by = ["Energy", "Area", "Ig"], ascending = [True, False, False])
    joined_df["Ig"] = 100 * joined_df["Ig"]
    joined_df["Ig"] = 100 * joined_df["Ig"].rename("Ig [%]")

    # joined_df.to_csv(f"with_Ig/{A}{element}.csv", index = False)
    return joined_df


def add_Ig(df, ig, ig_thr = 1.0):
    # UGLY, never for-loop in pandas
    added_df = df
    for row in df.iterrows():
        # WTF some pandas iterrows bizzare
        row = row[1]

        suitable_lines = ig[(row["Energy"] - row["FWHM"] < ig["E_tab"] + ig["sigm_E"]) & (row["Energy"] + row["FWHM"] > ig["E_tab"] - ig["sigm_E"]) & ig["Ig"] > 0.01*ig_thr]
        suitable_lines["Energy"] = float(row["Energy"])
        added_df = added_df.append(suitable_lines)

    print("Ig added")
    return added_df

if __name__ == "__main__":
    folder = argv[1]
    append_Igamma_dir(folder)