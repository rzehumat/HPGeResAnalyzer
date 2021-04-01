import rptParser
import getIg
import getEpsilon
import glob
import os
import json
import re

import pandas as pd
import uncertainties as uc

from collections import defaultdict
from pathlib import Path
from math import isnan
from addOrigin import addOrigin
from countRR import countRR
from getIg import permute_columns

OUTPUT_DIR = "out"


def add_si(num):
    return "\num{{{}}}".format(num)


def siunitx_mhchem(df):
    for col in ["Area", "RR_fiss_prod", "fiss_yield",
                "Half-life [s]"]:
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
    if isinstance(num, float):
        return num
    elif isnan(num.s):
        return num.n
    else:
        return '{:.1uS}'.format(num)


def polish_res(df):
    # UGLY, rewrite using for loops
    df["Energy"] = df[["Energy", "FWHM"]].apply(add_unc_en, axis=1)
    df["E_tab"] = df[["E_tab", "sigm_E"]].apply(add_unc, axis=1)
    df["fiss_yield"] = df[["fiss_yield", "sigm_fiss_yield"]].apply(
        add_unc, axis=1)
    df["old_fiss_yield"] = df["fiss_yield"].copy()
    df["Ig [%]"] = (100 * df["Ig"]).apply(unc_to_bracket)
    df["old_RR_fiss_prod"] = df["RR_fiss_prod"].copy()

    for col in ["RR", "Area", "RR_fiss_prod", "Real Time",
                "Live Time", "Half-life [s]"]:
        df[col] = df[col].apply(unc_to_bracket)

    return df


def drop_but_prodmode(df):
    sigm_cols = [x for x in df.columns.to_list() if "sigm_" in x]
    peak_analysis_report = [
        x for x in df.columns.to_list() if "Peak Analysis Report" in x]

    auxilliary = ["eps"]

    drop_cols = (sigm_cols
                 + ["FWHM", "%err", "Ig", "Dead Time",
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
    return df.drop(columns=drop_cols)


print("Available modes (default 0):")
print("0 ... 'Process dir' with terminal inputs"
      "(not recommended due to no error-safety)")
print("1 ... 'Process dir' with inputs from 'config.json' file (recommended)")

ig_all_df = pd.read_parquet("aux_data/ig_all.pq")
info_df = pd.read_parquet("aux_data/info_all.pq")
yield_df = pd.read_csv("aux_data/fissionYield_238U.csv", index_col=0)
mu_df = pd.read_csv("aux_data/mu_92.csv", index_col=0)
CALIBRATION_PATH = "aux_data/epsilons.csv"
eps_df = pd.read_csv(CALIBRATION_PATH, index_col=0)

print("Select mode:")
mode = input("[0]/[1]/[2] (default 1)") or "1"

if mode == "0":
    raw_dir = input("Enter relative path to directory"
                    "(press enter to keep default"
                    " 'raw_reports': ") or "raw_reports"
    if not os.path.isdir(raw_dir):
        raise Exception(f"Folder {raw_dir} not found."
                        "Consider creating it and moving RPT files there.")
    else:
        Path(OUTPUT_DIR).mkdir(parents=True, exist_ok=True)
        keys = ["A", "element", "detector_geometry", "foil_material_rho",
                "foil_thickness", "foil_mass",
                "foil_material_molar_mass", "irradiation_time",
                "irradiation_start"]
        units = ["mass number", "chemical abbr", "3/30/80/120/250",
                 "g.cm^-3", "cm", "g", "g.mol^-1", "s",
                 "dd.mm.yyyy hh:mm:ss"]

        kwargs = defaultdict(list)

        file_names = []
        for file_path in glob.iglob(f"{raw_dir}/*.RPT"):
            # UGLY, rewrite using the "unpacking operator *"
            # UGLY, reuse the loop-over-RPT-files
            file_names.append(file_path)
            try:
                file_vars = getEpsilon.file_name_parse(file_path)
                for i in range(3):
                    kwargs[keys[i]].append(file_vars[i])
            except ValueError:
                print("Info not detected from file name'"
                      f"{file_path.split('/')[-1]}'.")
                print("If A, element are known, add it manually."
                      "If not, leave A and element blank.")
                for i in range(3):
                    user_input = input(f"{keys[i]} [{units[i]}] = ")
                    kwargs[keys[i]].append(str(user_input))

            for i in range(3, len(keys)):
                if i < len(keys) - 1:
                    inp = input(f"{keys[i]} [{units[i]}]"
                                " (expected format '̈́0.23(1)') = ")
                    inp = uc.ufloat_fromstr(inp)
                    kwargs[keys[i]].append(inp)
                else:
                    kwargs[keys[i]].append(input(f"{keys[i]} [{units[i]}] = "))

        kwargs_df = pd.DataFrame.from_dict(kwargs)
        for j in range(len(file_names)):

            # Move to the dict above!!!
            print("Allowed time units [ns, us, ms, s, m, h, d, w, y]")
            print("Expected format '2.45 y'")
            hl_upper_bound = input(
                "Enter maximal Half-life to include (default 2 y) = ") or "2 y"
            hl_upper_bound = pd.to_timedelta(hl_upper_bound).total_seconds()

            ig_lower_bound = (input(
                "Enter minimal Ig to include [in %] (default 0.01 %) = ")
                 or "0.01")
            ig_lower_bound = float(ig_lower_bound)

            file_name = file_names[j]
            kwargs = kwargs_df.loc[j].to_dict()

            raw_df = rptParser.parse_one_RPT(file_name)
            df_ig = getIg.append_Igamma(raw_df, kwargs["A"],
                                        kwargs["element"], ig_all_df)
            df_ig_eps = getEpsilon.add_epsilon_file(
                df_ig, kwargs["detector_geometry"], eps_df)

            file_name = file_name.split("/")[-1].split(".")[-2]

            df_ig_eps_orig = addOrigin(df_ig_eps, info_df, yield_df)
            df_ig_eps_orig = countRR(df_ig_eps_orig, mu_df,
                                     kwargs["foil_material_rho"],
                                     kwargs["foil_thickness"],
                                     kwargs["foil_mass"],
                                     kwargs["foil_material_molar_mass"],
                                     kwargs["irradiation_time"],
                                     kwargs["irradiation_start"])
            # Half-life threshold should be
            df_ig_eps_orig = df_ig_eps_orig[
                (df_ig_eps_orig["Half-life [s]"] < hl_upper_bound)
                | (df_ig_eps_orig["FWHM"] > 0)]
            df_ig_eps_orig = df_ig_eps_orig[
                (df_ig_eps_orig["Ig [%]"] >= ig_lower_bound)
                | (df_ig_eps_orig["FWHM"] > 0)]
            prod_cols = [
                x for x in df_ig_eps_orig.columns.to_list() if "Prod_mode" in x]
            fff = ["Energy", "E_tab", "Ig [%]", "Area",
                   "Isotope", "RR", "RR_fiss_prod"]
            cols = fff + prod_cols
            df_ig_eps_orig = permute_columns(df_ig_eps_orig, cols)

            df_ig_eps_orig.to_csv(f"{OUTPUT_DIR}/{file_name}.csv", index=False)
            df_ig_eps_orig[
                (df_ig_eps_orig["FWHM"] > 0)  # to determine the original lines
                | (df_ig_eps_orig["Prod_mode_Fission product"])
                | (df_ig_eps_orig["fiss_yield"] > 0)
                ].to_csv(f"{OUTPUT_DIR}/{file_name}_fissile_products.csv",
                         index=False)
            df_ig_eps_orig[
                (df_ig_eps_orig["FWHM"] > 0)  # to determine the original lines
                | (df_ig_eps_orig["Prod_mode_Fast neutron activation"])
                | (df_ig_eps_orig["Prod_mode_Thermal neutron activation"])
                ].to_csv(f"{OUTPUT_DIR}/{file_name}_activation.csv",
                         index=False)

elif mode == "1":
    json_config = open("config.json", "r")
    config = json.load(json_config)
    # UGLY, reuse the loop below
    for file_path in glob.iglob(f"{config['raw_dir']}/*.RPT"):
        raw_df = rptParser.parse_one_RPT(file_path)
        jsn_config = open("config.json", "r")
        cfg = json.load(jsn_config)
        kwargs = get_kwargs(cfg, file_path)

        df = getIg.append_Igamma(raw_df, ig_all_df, **kwargs)
        print(f"file_path is {file_path}")
        df = getEpsilon.add_epsilon_file(df, eps_df, **kwargs)

        file_name = file_path.split("/")[-1].split(".")[-2]

        df = addOrigin(df, info_df, yield_df)
        df = countRR(df, mu_df, **kwargs)

        df.to_csv(f"{OUTPUT_DIR}/{file_name}_raw.csv", index=False)

        # format text output to human- and LaTeX-readable
        polished_df = polish_res(df)

        # drop all unnecessary columns but prodmode
        dropped_df = drop_but_prodmode(polished_df)

        # restrict hl and Ig interval
        # hl_upper_bound = pd.to_timedelta(kwargs['hl_upper_bound']).total_seconds()
        hl_lower_bound = pd.to_timedelta(
            kwargs['hl_lower_bound']).total_seconds()

        dropped_df.to_csv(
            f"{OUTPUT_DIR}/{file_name}_polished.csv", index=False)

        fiss_df = dropped_df[
            (dropped_df["E_tab"].isna())  # to determine the original lines
            | (dropped_df["old_RR_fiss_prod"] > 0)
            ]
        prod_cols = [x for x in fiss_df.columns.to_list() if "Prod_mode" in x]
        fiss_df = fiss_df.drop(columns=prod_cols)

        first_cols = ["Energy", "E_tab", "Ig [%]", "Area", "Isotope",
                      "RR", "RR_fiss_prod", "fiss_yield", "Half-life [s]"]
        fiss_df = permute_columns(fiss_df, first_cols)

        fiss_df = fiss_df[
            (fiss_df["old_RR_fiss_prod"] > 1e-17)
            & (fiss_df["old_RR_fiss_prod"] < 1e-14)]

        fiss_df.to_csv(f"{OUTPUT_DIR}/{file_name}_fissile_products.csv",
                       index=False)
        fiss_df = siunitx_mhchem(fiss_df)
        fiss_df["Isotope"] = fiss_df["Isotope"].apply(to_mhchem)
        fiss_df_tex = fiss_df.to_latex(index=False,
                                       columns=["Energy", "E_tab", "Ig [%]",
                                                "Area", "Isotope",
                                                "RR_fiss_prod", "fiss_yield",
                                                "Half-life [s]"],
                                       buf=f"{OUTPUT_DIR}/{file_name}_fissile_products.tex",
                                       na_rep="", position="h",
                                       caption=(f"{file_name}", ""),
                                       label=f"{file_name}",
                                       escape=False,
                                       longtable=True,
                                       column_format="SSSSlSSS")

        # df = df[
        #         (df["Half-life [s]"] >= hl_lower_bound)
        #         & (df["Ig [%]"] >= kwargs["ig_lower_bound"])
        #         & (df["Ig [%]"] <= kwargs["ig_upper_bound"])
        #         | (df["FWHM"] > 0)]

        # prod_cols = [x for x in df.columns.to_list() if "Prod_mode" in x]
        # fff = ["Energy", "E_tab", "Ig [%]", "Area",
        #        "Isotope", "RR", "RR_fiss_prod"]
        # cols = fff + prod_cols
        # df = permute_columns(df, cols)
        # df.to_csv(f"{OUTPUT_DIR}/{file_name}.csv", index=False)
        # df[
        #     (df["FWHM"] > 0)  # to determine the original lines
        #     # | (df["Prod_mode_Fission product"])
        #     | (df["fiss_yield"] > 0)
        #     ].to_csv(f"{OUTPUT_DIR}/{file_name}_fissile_products.csv",
        #              index=False)
        # df[
        #     (df["FWHM"] > 0)  # to determine the original lines
        #     | (df["Prod_mode_Fast neutron activation"])
        #     | (df["Prod_mode_Thermal neutron activation"])
        #     ].to_csv(f"{OUTPUT_DIR}/{file_name}_activation.csv",
        #              index=False)
