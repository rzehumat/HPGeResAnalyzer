import glob
import os
import pandas as pd
from pathlib import Path
import re
import dateutil.parser as dparser


def parse_header(lines):
    header = {}
    datetime_columns = []

    for line in lines:
        if not line.isspace():
            if ":" in line:
                new_key, new_val = [x.strip() for x in line.split(':',maxsplit=1)]
                if new_val == "":
                    continue
                if re.match("\d{1,2}\.\d{1,2}\.\d{2,4}[ \d{1,2}\:\d{1,2}\:\d{1,2}]*", new_val):
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


def parse_RPT(rpt_file):
    my_file = open(rpt_file, 'r')
    lines = my_file.readlines()

    lookup = "Bkgnd"
    for line in lines:
        if lookup in line:
            split_num = lines.index(line)
            break

    header_lines = lines[0:split_num]
    data_lines = lines[split_num:-1]

    header, datetime_columns = parse_header(header_lines)
    data = parse_data(data_lines)

    meta_df = pd.DataFrame(header, index=[0])
    rpt_df = pd.DataFrame(data)
    res = pd.concat([rpt_df, meta_df], axis=1)
    res = res.fillna(method='ffill')
    
    return res


if not os.path.isdir("raw_reports"):
    raise Exception("Folder 'raw_reports' not found. Consider creating it and moving RPT files there.")

Path("./parsed_reports").mkdir(parents=True, exist_ok=True)

for rpt_file in glob.iglob("./raw_reports/*.RPT"):
    res_df = parse_RPT(rpt_file)

    name = (rpt_file.split('.')[-2]).split('/')[-1]
    res_df.to_csv(f"parsed_reports/{name}.csv")