def addOrigin(df, info, yield_df):
    cols = info.columns.to_list()

    df[cols] = info[cols]
    df["fiss_yield"] = yield_df["Yield"]
    df["sigm_fiss_yield"] = yield_df["Error"]

    return df
