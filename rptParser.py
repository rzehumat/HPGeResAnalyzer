import glob
import os
import re
import datetime
import pandas as pd
import dateutil.parser as dparser
from sys import argv
from pathlib import Path



def parse_header(lines):
    header = {}
    datetime_columns = []

    for line in lines:
        if not line.isspace():
            if ":" in line:
                new_key, new_val = [x.strip() for x in line.split(':',maxsplit=1)]
                if new_val == "":
                    continue
                if re.match(r"\d{1,2}\.\d{1,2}\.\d{2,4}[ \d{1,2}\:\d{1,2}\:\d{1,2}]*", new_val):
                    new_val = dparser.parse(new_val, fuzzy=True)
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
        if not "Bkgnd" in line:
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
    
    return res

def parse_RPT(folder):
    if not os.path.isdir(f"./{folder}"):
        raise Exception(f"Folder {folder} not found. Consider creating it and moving RPT files there.")

    Path("./parsed_reports").mkdir(parents=True, exist_ok=True)

    for rpt_file in glob.iglob(f"./{folder}/*.RPT"):
        res_df = parse_one_RPT(rpt_file)

        name = (rpt_file.split('.')[-2]).split('/')[-1]

        cols = res_df.columns.tolist()
        for col in DATETIME_COLUMNS:
            print(col)
            if col in cols:
                print("AAA")
                res_df[col] = pd.to_pydatetime(res_df[col])

        res_df.to_csv(f"parsed_reports/{name}.csv")

if __name__ == "__main__":
    folder = argv[1]

    DATETIME_COLUMNS = ["Report Generated On"]

    parse_RPT(folder)