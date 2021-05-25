from os import error
import rptParser
import getIg
import getEpsilon
import glob
# import os
import json
import re

import pandas as pd
import numpy as np
import uncertainties as uc

# from collections import defaultdict
# from pathlib import Path
from math import isnan
from addOrigin import addOrigin
from countRR import countRR
from getIg import permute_columns
from uncertainties import unumpy

OUTPUT_DIR = "out"


def average_rr(rr_uc_df, kind):
    rr = pd.DataFrame(unumpy.nominal_values(rr_uc_df))
    std = pd.DataFrame(unumpy.std_devs(rr_uc_df))
    std[0] = std[0].fillna(0.01 * rr[0])
    std = std[~rr.isna()].dropna()
    rr = rr.dropna()
    weights = 1 / (std * std)
    RR = (weights * rr).sum()/weights.sum()
    std_RR = np.sqrt(1/weights.sum())
    result = add_si(add_unc([RR.loc[0], std_RR.loc[0]]))

    out_file = open(f"{OUTPUT_DIR}/{file_name}_{kind}_result.txt", "w")
    out_file.write(result)
    out_file.close()


def add_si(num):
    return "\\num{{{}}}".format(num)


def siunitx_mhchem(df):
    for col in ["Area", "RR", "Half-life [s]"]:
        df[col] = df[col].apply(add_si)
    return df


def to_mhchem(isotope_str):
    element = re.findall(r'[A-Z]{1}.*', isotope_str)[0]
    A = isotope_str.replace(element, "")
    return "\ce{{^{{{}}}{}}}".format(A, element)


def add_unc(x):
    return '{:.1uS}'.format(uc.ufloat(x[0], x[1]))


def get_kwargs(config, file_path):
    my_dict = config
    for key in list(config.keys()):
        if key not in file_path:
            my_dict.pop(key, None)
    df = pd.json_normalize(my_dict, sep='.')
    mydict = df.to_dict(orient='records')[0]
    # Ugly
    for key in tuple(mydict.keys()):
        mydict[key.split('.')[-1]] = mydict.pop(key)
    return mydict


def add_unc_en(x):
    return '{:.1uS}'.format(uc.ufloat(x[0], 0.5*x[1]))


def unc_to_bracket(num):
    if isinstance(num, float) or isinstance(num, int):
        return num
    elif isnan(num.s):
        num = uc.ufloat(num.n, 0.01*num.n)
        return '{:.1uS}'.format(num)
    else:
        return '{:.1uS}'.format(num)


def polish_res(df):
    # UGLY, rewrite using for loops
    df["old_Energy"] = df["Energy"].copy()
    df["old_Area"] = df["Area"].copy()
    df["Energy"] = df[["Energy", "FWHM"]].apply(add_unc_en, axis=1)
    df["E_tab"] = df[["E_tab", "sigm_E"]].apply(add_unc, axis=1)
    df["fiss_yield"] = df[["fiss_yield", "sigm_fiss_yield"]].apply(
        add_unc, axis=1)
    df["old_fiss_yield"] = df["fiss_yield"].copy()
    
    if kwargs["RR"]:
        df["old_RR"] = df["RR"].copy()
    df["Ig [%]"] = (100 * df["Ig"]).apply(unc_to_bracket)
    if kwargs["RR"]:
        df["old_RR_fiss_prod"] = df["RR_fiss_prod"].copy()

    cols = ["Area", "Real Time", "Live Time", "Half-life [s]"]
    if kwargs["RR"]:
        cols += ["RR", "RR_fiss_prod"]
    for col in cols:
        df[col] = df[col].apply(unc_to_bracket)

    return df


def drop_but_prodmode(df):
    sigm_cols = [x for x in df.columns.to_list() if "sigm_" in x]
    peak_analysis_report = [
        x for x in df.columns.to_list() if "Peak Analysis Report" in x]

    auxilliary = ["eps"]

    drop_cols = (sigm_cols
                 + ["FWHM", "%err", "Dead Time",
                    "Peak Locate Threshold", "Sample Geometry",
                    "Use Fixed FWHM", "Pk", "IT", "Bkgnd", "Channel",
                    "Left", "PW", "Cts/Sec", "Fit", "Filename",
                    "Sample Identification", "Sample Type",
                    "Peak Locate Range (in channels)",
                    "Peak Area Range (in channels)",
                    "Identification Energy Tolerance", "Sample Size",
                    "Efficiency ID", "Peak Analysis From Channel",
                    "Peak Search Sensitivity", "Max Iterations",
                    "Peak Fit Engine Name", "Stable", "Update",
                    "Literature cut-off date", "Author(s)", "E(level)",
                    "Sn(keV)", "References since cut-off", "Sp(keV)",
                    "Abundance", "Jp", "ENSDF citation",
                    "Energy Calibration Used Done On",
                    "Efficiency Calibration Used Done On"]
                 + auxilliary
                 + peak_analysis_report)
    return df.drop(columns=drop_cols, errors="ignore")


print("Available modes (default 1):")
print("1 ... 'Process dir' with inputs from 'config.json' file (recommended)")
print("2 ... Test mode")

ig_all_df = pd.read_parquet("aux_data/ig_all.pq")
info_df = pd.read_parquet("aux_data/info_all.pq")
yield_df = pd.read_csv("aux_data/fissionYield_238U.csv", index_col=0)
mu_df = pd.read_csv("aux_data/mu_92.csv", index_col=0)
CALIBRATION_PATH = "aux_data/epsilons.csv"
eps_df = pd.read_csv(CALIBRATION_PATH, index_col=0)

print("Select mode:")
mode = input("[1]/[2] (default 1)") or "1"

# if mode == "1":
if mode == "1":
    conf_file_name = "config.json"
elif mode == "2":
    conf_file_name = "test/config_test.json"

# json_config = open("config.json", "r")
json_config = open(conf_file_name, "r")
config = json.load(json_config)
# UGLY, reuse the loop below
for file_path in glob.iglob(f"{config['raw_dir']}/*.RPT"):
    raw_df = rptParser.parse_one_RPT(file_path)
    # jsn_config = open("config.json", "r")
    jsn_config = open(conf_file_name, "r")
    cfg = json.load(jsn_config)
    kwargs = get_kwargs(cfg, file_path)

    df = getIg.append_Igamma(raw_df, ig_all_df, **kwargs)
    print(f"file_path is {file_path}")
    df = getEpsilon.add_epsilon_file(df, eps_df, **kwargs)

    file_name = file_path.split("/")[-1].split(".")[-2]

    # if kwargs["RR"]:
    #     print("AAA")
    # else:
    #     print("BBB")
    # input()

    df = addOrigin(df, info_df, yield_df)
    # Sometimes we are not interested in reaction rates
    if kwargs["RR"]:
        df = countRR(df, mu_df, **kwargs)

    df.to_csv(f"{OUTPUT_DIR}/{file_name}_raw.csv", index=False)

    hl_lower_bound = pd.to_timedelta(
        kwargs['hl_lower_bound']).total_seconds()

    df = df[df["Half-life [s]"] > hl_lower_bound]

    # format text output to human- and LaTeX-readable
    polished_df = polish_res(df)

    # drop all unnecessary columns but prodmode
    dropped_df = drop_but_prodmode(polished_df)

    activation_df = dropped_df[
            (dropped_df["Isotope"] == "238U")
            | (df["Isotope"] == "239U")
            | (df["Isotope"] == "239Np")
            | (df["Isotope"] == "239Pu")]
    # activation_df = activation_df[activation_df["Ig"] > 0.00001]

    prod_cols = [x for x in activation_df.columns.to_list()
                 if "Prod_mode" in x]
    activation_df = activation_df.drop(columns=prod_cols)

    first_cols = ["Energy", "E_tab", "Ig [%]", "Area", "Isotope",
                  "fiss_yield", "Half-life [s]"]
    if kwargs["RR"]:
        first_cols += ["RR", "RR_fiss_prod"]
    activation_df = permute_columns(activation_df, first_cols)

    activation_df.to_csv(f"{OUTPUT_DIR}/{file_name}_activation.csv",
                         index=False)
    activation_df = siunitx_mhchem(activation_df)
    activation_df["Isotope"] = activation_df["Isotope"].apply(to_mhchem)

    dropped_df.to_csv(
        f"{OUTPUT_DIR}/{file_name}_polished.csv", index=False)

    fiss_df = dropped_df[
        (dropped_df["E_tab"].isna())  # to determine the original lines
        | (dropped_df["old_RR_fiss_prod"] > 0)
        ]
    fiss_df = fiss_df[fiss_df["Ig"] > 0.01]
    prod_cols = [x for x in fiss_df.columns.to_list() if "Prod_mode" in x]
    fiss_df = fiss_df.drop(columns=prod_cols)

    first_cols = ["Energy", "E_tab", "Ig [%]", "Area", "Isotope",
                  "RR", "RR_fiss_prod", "fiss_yield", "Half-life [s]"]
    fiss_df = permute_columns(fiss_df, first_cols)

    fiss_df = fiss_df[
        (fiss_df["old_RR_fiss_prod"] > 1e-18)
        & (fiss_df["old_RR_fiss_prod"] < 1e-14)]

    dh = pd.DataFrame()

    for energy in fiss_df["old_Energy"].unique():
        dg = fiss_df[fiss_df["old_Energy"] == energy]

        dg["rel_fiss_RR"] = 1e+16 * dg["old_RR_fiss_prod"]
        dg["prod_exc"] = dg["rel_fiss_RR"].product() / dg["rel_fiss_RR"]
        dg["fraction"] = dg["prod_exc"] / dg["prod_exc"].sum()

        dg["mod_fiss_RR"] = dg["fraction"] * dg["old_RR_fiss_prod"]

        dg["mod_Area"] = dg["fraction"] * dg["old_Area"]
        dg["mod_Area"] = dg["mod_Area"].apply(unc_to_bracket)

        dh = dh.append(dg)

    fiss_df["mod_fiss_RR"] = dh["mod_fiss_RR"]
    average_rr(fiss_df["mod_fiss_RR"], "fiss")
    dg["mod_fiss_RR"] = dg["mod_fiss_RR"].apply(unc_to_bracket)
    fiss_df["mod_Area"] = dh["mod_Area"]

    fiss_df.to_csv(f"{OUTPUT_DIR}/{file_name}_fissile_products.csv",
                   index=False)
    fiss_df = siunitx_mhchem(fiss_df)
    fiss_df["Isotope"] = fiss_df["Isotope"].apply(to_mhchem)

    activation_df = activation_df[
        ~activation_df["old_Energy"].isin(fiss_df["old_Energy"])]
    average_rr(activation_df["old_RR"], "activation")
    activation_df["Ig [%]"] = activation_df["Ig [%]"].apply(add_si)
    activation_df_tex = activation_df.to_latex(
        index=False,
        columns=["Energy", "E_tab", "Ig [%]", "Area", "Isotope",
                 "RR", "Half-life [s]"],
        buf=f"{OUTPUT_DIR}/{file_name}_activation.tex",
        header=["{$E\_{mer}$ [\si{keV}]}", "{$E\_{tab}$ [\si{keV}]}",
                "$I_{\gamma}$ [\%]", "$S\_{peak}$ [\si{keV}]",
                "Izotop", "$RR_{(n,g)}$ [\si{cm^{-3} s^{-1}}]",
                "$T_{1/2}$ [\si{s}]"],
        na_rep="", position="h",
        caption=(f"{file_name}", ""),
        label=f"{file_name[:5]}-rc",
        escape=False,
        longtable=True,
        column_format="SScclcc")

    fiss_df["mod_fiss_RR"] = fiss_df["mod_fiss_RR"].apply(unc_to_bracket)
    fiss_df["mod_fiss_RR"] = fiss_df["mod_fiss_RR"].apply(add_si)
    fiss_df["mod_Area"] = fiss_df["mod_Area"].apply(add_si)
    fiss_df_tex = fiss_df.to_latex(
        index=False,
        columns=["Energy", "E_tab", "Ig [%]", "Area", "mod_Area",
                 "Isotope", "mod_fiss_RR", "fiss_yield", "Half-life [s]"],
        header=["{$E\_{mer}$ [\si{keV}]}",
                "{$E\_{tab}$ [\si{keV}]}",
                "$I_{\gamma}$ [\%]",
                "$S\_{peak}$ [\si{keV}]",
                "$S\_{mod}$ [\si{keV}]",
                "Izotop",
                "\\footnotesize{$\\frac{RR_{(n,f)}}{[\si{cm^{-3} s^{-1}}]}$}",
                "$Y_f$",
                "$T_{1/2}$ [\si{s}]"],
        buf=f"{OUTPUT_DIR}/{file_name}_fissile_products.tex",
        na_rep="", position="h",
        caption=(f"{file_name}", ""),
        label=f"{file_name[:5]}-f",
        escape=False,
        longtable=True,
        column_format="Scccrlrrr")
