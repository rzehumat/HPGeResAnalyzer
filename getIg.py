import pandas as pd
import os
import glob

from pathlib import Path
from sys import argv
from getEpsilon import file_name_parse


def append_Igamma_dir(parsed_dir):
    if not os.path.isdir(parsed_dir):
        raise Exception(f"Directory '{parsed_dir}' not found."
                        "Consider creating it and moving parsed files there.")

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


def append_Igamma(parsed_df, ig_all_df, **kwargs):
    print(kwargs)
    if "A" in kwargs.keys():
        ig_df = ig_all_df.loc[
            [f"{kwargs['A']}{kwargs['element']}",
             f"{int(kwargs['A'])+1}{kwargs['element']}",
             f"{int(kwargs['A'])+2}{kwargs['element']}"]]
    else:
        ig_df = ig_all_df

    joined_df = add_Ig(parsed_df, ig_df)

    Path("./with_Ig").mkdir(parents=True, exist_ok=True)
    joined_df = joined_df.sort_values(by=["Energy", "Pk", "Ig"],
                                      ascending=[True, False, False])
    joined_df["Ig [%]"] = 100 * joined_df["Ig"].rename("Ig [%]")
    joined_df = permute_columns(joined_df, ["Energy", "E_tab", "Ig [%]",
                                            "Area", "Isotope", "sigm_E",
                                            "sigm_Ig", "FWHM", "%err",
                                            "Live Time", "Real Time",
                                            "Dead Time (rel)"])

    return joined_df


def add_Ig(df, ig, ig_thr=1.0):
    # UGLY, never for-loop in pandas
    added_df = df
    print(df.shape[0])
    for row in range(df.shape[0]):
        condition = (
            (df.loc[row, "Energy"] - df.loc[row, "FWHM"] < ig["E_tab"] + ig["sigm_E"])
            & (df.loc[row, "Energy"] + df.loc[row, "FWHM"] > ig["E_tab"] - ig["sigm_E"])
            & ig["Ig"] > 0.01*ig_thr)

        suitable_lines = ig[condition]
        suitable_lines["Isotope"] = ig[condition].index
        suitable_lines.loc[:, ("Energy")] = df.loc[row, "Energy"]
        suitable_lines.loc[:, ("Area")] = df.loc[row, "Area"]
        suitable_lines.loc[:, ("FWHM")] = df.loc[row, "FWHM"]
        suitable_lines.loc[:, ("%err")] = df.loc[row, "%err"]

        added_df = added_df.append(suitable_lines)

    print("Ig added")
    return added_df


if __name__ == "__main__":
    folder = argv[1]
    append_Igamma_dir(folder)
