import pandas as pd


def addOrigin(df, info):
    cols_drop = ["Update", "Sn(keV)", "sigm_Sn(keV)", "E(level)", "sigm_E(level)", "Literature cut-off date", "Author(s)", "Sp(keV)", "sigm_Sp(keV)", "References since cut-off", "Jp", "ENSDF citation"]
    cols = [item for item in info.columns.to_list() if item not in cols_drop]

    df[cols] = info[cols]

    return df