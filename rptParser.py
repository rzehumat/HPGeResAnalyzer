import glob
import os
import re

import pandas as pd
import dateutil.parser as dparser

from sys import argv
from pathlib import Path
from getIg import permute_columns
from uncertainties import ufloat_fromstr


def parse_header(lines):
    header = {}
    datetime_columns = []

    for line in lines:
        if not line.isspace():
            if ":" in line:
                new_key, new_val = [
                    x.strip() for x in line.split(':', maxsplit=1)
                    ]
                if new_val == "":
                    continue
                if re.match(
                    r"\d{1,2}\.\d{1,2}\.\d{2,4}[ \d{1,2}\:\d{1,2}\:\d{1,2}]*",
                   new_val):
                    new_val = dparser.parse(new_val, fuzzy=True, dayfirst=True)
                    datetime_columns.append(new_key)

                header[new_key] = new_val

    return header, datetime_columns


def parse_data(lines):
    # UGLY, could be one-liner, or?
    head = lines[0].split()
    params_no = len(head)

    data = dict.fromkeys(head)

    # UGLY, I believe there shoould be a way of doing it without for loop
    for key in data:
        data[key] = []

    for line in lines[1:]:
        if "Bkgnd" not in line:
            if len(line.split()) == params_no:
                # UGLY, do it in one line
                for i in range(params_no):
                    data[head[i]].append(line.split()[i])

    return data


def parse_one_RPT(rpt_file):
    my_file = open(rpt_file, 'r')
    lines = my_file.readlines()

    lookup = "Bkgnd"
    for line in lines:
        if lookup in line:
            split_num = lines.index(line)
            break

    header_lines = lines[0:split_num]
    data_lines = lines[split_num:-1]

    header, _ = parse_header(header_lines)
    data = parse_data(data_lines)

    meta_df = pd.DataFrame(header, index=[0])
    rpt_df = pd.DataFrame(data)
    res = pd.concat([rpt_df, meta_df], axis=1)
    res = res.fillna(method='ffill')

    res = polish_dtypes(res)
    # res = permute_columns(res, [
    #                             "Pk", "Energy", "Area", "Live Time",
    #                             "Real Time", "Dead Time (rel)",
    #                             "Acquisition Started", "Bkgnd",
    #                             "FWHM", "Channel", "Cts/Sec", "%err"
    #                             ])
    return res


def add_uncert(t):
    return ufloat_fromstr(str(t) + "0(5)")


def polish_dtypes(df):
    DATETIME_COLUMN = "Report Generated On"
    UNSIGNED_COLUMNS = ["Area", "Bkgnd", "Left", "PW", "Sample Identification"]
    INT_COLUMNS = ["IT", "Sample Type"]
    FLOAT_COLUMNS = [
        "Energy", "FWHM", "Channel", "Cts/Sec", "%err", "Fit",
        "Peak Locate Threshold"
        ]
    # TO_DROP = [
    #     "Sample Geometry", "Peak Locate Range (in channels)", "Sample Size",
    #     "Dead Time",
    #     "Peak Analysis Report                    26.11.2020  5",
    #     "Peak Analysis From Channel", "Peak Search Sensitivity",
    #     "Max Iterations",
    #     "Use Fixed FWHM", "Peak Fit Engine Name", "Left", "PW",
    #     "Fit", "Filename",
    #     "Sample Identification", "Sample Type", "Peak Locate Threshold",
    #     "Peak Area Range (in channels)", "Efficiency ID"]
    cols = df.columns.tolist()
    if DATETIME_COLUMN in cols:
        df[DATETIME_COLUMN] = pd.to_datetime(df[DATETIME_COLUMN])
    for col in UNSIGNED_COLUMNS:
        df[col] = pd.to_numeric(df[col], downcast="unsigned")
    for col in INT_COLUMNS:
        df[col] = pd.to_numeric(df[col], downcast="integer")
    for col in FLOAT_COLUMNS:
        df[col] = pd.to_numeric(df[col], downcast="float")

    col_name = "Identification Energy Tolerance"
    if col_name in cols:
        unit = str(df[col_name][0]).split()[-1]
        df[col_name] = (df[col_name].str.split(" ").str[0]).astype(float)
        df[col_name] = df[col_name].rename(f"{col_name} [{unit}]")
    for col in ["Real Time", "Live Time"]:
        if col in cols:
            unit = str(df[col][0]).split()[-1]
            df[col] = df[col].rename(f"{col} [{unit}]")
            df[col] = df[col].str.split(" ").str[0].apply(add_uncert)

    if "Dead Time" in cols:
        df["Dead Time [%]"] = (
            df["Dead Time"].str.split(" ").str[0]
            ).astype(float)

    # df = df.drop(columns=TO_DROP, errors="ignore")
    return df


def parse_RPT(folder):
    if not os.path.isdir(f"./{folder}"):
        raise Exception(f"Folder {folder} not found."
                        "Consider creating it and moving RPT files there.")

    Path("./parsed_reports").mkdir(parents=True, exist_ok=True)

    for rpt_file in glob.iglob(f"./{folder}/*.RPT"):
        res_df = parse_one_RPT(rpt_file)

        name = (rpt_file.split('.')[-2]).split('/')[-1]

        res_df.to_csv(f"parsed_reports/{name}.csv")


if __name__ == "__main__":
    folder = argv[1]
    DATETIME_COLUMNS = ["Report Generated On"]
    parse_RPT(folder)
