import pandas as pd
import numpy as np

import uncertainties.unumpy as unp
from uncertainties import ufloat_fromstr
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


def parse_time_unc(time_str):
    num_str, unit = time_str.split(' ')
    num = ufloat_fromstr(num_str)
    return num * pd.to_timedelta(f"1 {unit}").total_seconds()


# def countRR(orig_df, mu_df, rho, d, mass, molar_mass,
#             t_irr_str, irr_start_str):
def countRR(orig_df, mu_df, **kwargs):
    """
    Assumption: One does not measure an isotope after 10 half-lives.
    Therefore we calculate delta_t and compare it with half-lives.
    If delta_t < 10HL => select df[Half-life] > 0.1 delta_t
    """
    irr_start = pd.to_datetime(kwargs["irradiation_start"], dayfirst=True)
    # irr_start = int(time.mktime(irr_start.timetuple()))

    t_irr = parse_time_unc(kwargs["irradiation_time"])
    # irr_end = irr_start + t_irr

    acq_started = pd.to_datetime(orig_df["Acquisition Started"][0],
                                 dayfirst=True)
    delta_t = (acq_started - irr_start).total_seconds()
    # delta_t = int(time.mktime(acq_started.timetuple())) - irr_end
    # delta_t = acq_started.total_seconds() - irr_start.total_seconds()

    df = orig_df[
        (orig_df["Half-life [s]"] > 0.1 * delta_t)
        & (orig_df["Half-life [s]"] < pd.to_timedelta("2 y").total_seconds())]
    lines_df = orig_df[orig_df["FWHM"] > 0]

    AVOGADRO = float(6.02214076e+23)
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

    df["RR"] = ((real_time / live_time)
                * k * lam * df["Area"].to_numpy()
                / (
                    N * (1-unp.exp(-lam*t_irr)) * unp.exp(-lam * delta_t)
                    * (1-unp.exp(-lam * real_time))
                    * df["eps"] * df["Ig"]))
    df["RR_fiss_prod"] = (2 / df["fiss_yield"]) * df["RR"]
    print("Reaction rates counted successfully.")
    df = df.append(lines_df)
    df = df.sort_values(by=["Energy", "FWHM", "Ig [%]"],
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
