import pandas as pd
import numpy as np
import time

# from uncertainties import ufloat
import uncertainties.unumpy as unp
from sys import argv


def get_mu_one(val, mu_df):
    if val in mu_df["photon_energy"]:
        return mu_df.loc[val, "total_attenuation_without_coherent_scatt"]
    lower = mu_df[mu_df["photon_energy"] < val].iloc[-1]
    upper = mu_df[mu_df["photon_energy"] > val].iloc[0]
    return (lower["total_attenuation_without_coherent_scatt"]
            + ((val-lower["photon_energy"])/((upper-lower)["photon_energy"]))
            * ((upper-lower)["total_attenuation_without_coherent_scatt"]))


def get_mu_col(energy_col, mu_df):
    return energy_col.apply(lambda x: get_mu_one(x, mu_df))


def countRR(orig_df, mu_df, rho, d, mass, molar_mass, t_irr, irr_start_str):
    """
    Assumption: One does not measure an isotope after 10 half-lives.
    Therefore we calculate delta_t and compare it with half-lives.
    If delta_t < 10HL => select df[Half-life] > 0.1 delta_t
    """
    irr_start = pd.to_datetime(irr_start_str, dayfirst=True)
    irr_start = int(time.mktime(irr_start.timetuple()))
    irr_end = irr_start + t_irr

    acq_started = pd.to_datetime(orig_df["Acquisition Started"][0])
    delta_t = int(time.mktime(acq_started.timetuple())) - irr_end

    df = orig_df[
        (orig_df["Half-life [s]"] > 0.1 * delta_t)
        | (orig_df["FWHM"] > 0)]

    AVOGADRO = float(6.02214076e+23)
    mu = get_mu_col(df["Energy"], mu_df)

    mu = mu.astype(np.float64)

    k = (mu * rho * d) / (1 - unp.exp(- mu * rho * d))
    lam = np.log(2) / df["Half-life [s]"]
    N = mass*AVOGADRO/molar_mass

    df["RR"] = ((df['Real Time'][0] / df["Live Time"][0])
                * k * lam * df["Area"].to_numpy()
                / (
                    N * (1-unp.exp(-lam*t_irr)) * unp.exp(-lam * delta_t)
                    * (1-unp.exp(-lam * df["Real Time"][0]))
                    * df["eps"] * df["Ig"]))

    print("Reaction rates counted successfully.")
    return df

    # def calc_uexp(lam_df, delta_t, N):
    #     length = lam_df.shape[0]
    #     lam_arr = lam_df.to_numpy()
    #     res = np.ones(length)
    #     for i in range(length):
    #         try:
    #             res[i] = unp.exp(lam_arr[i] * delta_t) / N
    #         except OverflowError:
    #             print(f"Overflow at i = {i} which is lam * delta = {lam_arr[i]*delta_t}")
    #             res[i] =
    #     return res


if __name__ == "__main__":
    path = argv[1]
    df = pd.read_csv(path)
    rho = 19.01
    d = 0.028
    mass = 0.2361
    molar_mass = 238.0507847
    mu_df = pd.read_csv("aux_data/mu_92.csv", index_col=0)
    t_irr = 30*60
    irr_start_str = "1.12.2020 10:31:00"
    countRR(df, mu_df, rho, d, mass, molar_mass, t_irr, irr_start_str)
    pd.to_csv(f"{path.split('.')[-1]}_wRR.csv")
