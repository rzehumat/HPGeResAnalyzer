import pandas as pd
import numpy as np
import uncertainties as uc

import uncertainties.unumpy as unp
from uncertainties import ufloat_fromstr
from sys import argv


def uncert_series(x):
    return uc.ufloat(x[0], x[1])


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


def parse_time_unc(time_str):
    num_str, unit = time_str.split(' ')
    num = ufloat_fromstr(num_str)
    return num * pd.to_timedelta(f"1 {unit}").total_seconds()


def countRR(orig_df, mu_df, **kwargs):
    """
    Assumption: One does not measure an isotope after 10 half-lives.
    Therefore we calculate delta_t and compare it with half-lives.
    If delta_t < 10HL => select df[Half-life] > 0.1 delta_t
    """
    AVOGADRO = float(6.02214076e+23)

    orig_df["Half-life [s]"] = orig_df[
        ["Half-life [s]", "sigm_Half-life [s]"]].apply(uncert_series, axis=1)
    orig_df["sigm_Area"] = 0.01 * orig_df["Area"] * orig_df["%err"]
    orig_df["Area"] = orig_df[
        ["Area", "sigm_Area"]].apply(uncert_series, axis=1)
    orig_df["Ig"] = orig_df[
        ["Ig", "sigm_Ig"]].apply(uncert_series, axis=1)

    irr_start = pd.to_datetime(kwargs["irradiation_start"], dayfirst=True)
    t_irr = parse_time_unc(kwargs["irradiation_time"])
    acq_started = pd.to_datetime(orig_df["Acquisition Started"][0],
                                 dayfirst=True)
    delta_t = (acq_started - irr_start).total_seconds() - t_irr
    T_LOW = 0.1 * delta_t
    T_HIGH = 6e+6

    df_low = orig_df[(orig_df["Half-life [s]"] < T_LOW)]
    df = orig_df[
        (orig_df["Half-life [s]"] > T_LOW)
        & (orig_df["Half-life [s]"] < T_HIGH)]
    lines_df = orig_df[orig_df["E_tab"].isna()]
    df_high = orig_df[(orig_df["Half-life [s]"] >= T_HIGH)]

    mu = get_mu_col(df["Energy"], mu_df)
    mu = mu.astype(np.float64)

    rho = ufloat_fromstr(kwargs['foil_material_rho'])
    d = ufloat_fromstr(kwargs['foil_thickness'])
    k = (mu * rho * d) / (1 - unp.exp(- mu * rho * d))

    lam = np.log(2) / df["Half-life [s]"]

    mass = ufloat_fromstr(kwargs['foil_mass'])
    molar_mass = ufloat_fromstr(kwargs['foil_material_molar_mass'])

    N = mass*AVOGADRO/molar_mass

    real_time = lines_df['Real Time'][0]
    live_time = lines_df['Live Time'][0]

    nom = ((real_time / live_time) * k * lam * df["Area"])
    denom = (N * (1-unp.exp(-lam*t_irr)) * unp.exp(-lam * delta_t)
             * (1-unp.exp(-lam * real_time))
             * df["eps"] * df["Ig"])

    df["RR"] = nom / denom

    df["RR_fiss_prod"] = (2 / df["fiss_yield"]) * df["RR"]
    print("Reaction rates counted successfully.")
    df = df.append(lines_df)
    # df = df.append(df_low)
    df = df.append(df_high)
    df = df.sort_values(by=["Energy", "Channel", "Ig [%]"],
                        ascending=[True, True, False])
    return df


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
