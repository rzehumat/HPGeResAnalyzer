def addOrigin(df, info, yield_df):
    # cols_drop = [
    #     "Update", "Sn(keV)", "sigm_Sn(keV)", "E(level)",
    #     "sigm_E(level)", "Literature cut-off date", "Author(s)",
    #     "Sp(keV)", "sigm_Sp(keV)", "References since cut-off",
    #     "Jp", "ENSDF citation"]
    # cols = [item for item in info.columns.to_list() if item not in cols_drop]
    cols = info.columns.to_list()

    df[cols] = info[cols]
    df["fiss_yield"] = yield_df["Yield"]
    df["sigm_fiss_yield"] = yield_df["Error"]

    return df
