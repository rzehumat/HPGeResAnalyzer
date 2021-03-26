import pandas as pd
import numpy as np

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


def countRR(df, mu_df, rho, d, mass, molar_mass, t_irr, irr_start_str):
    AVOGADRO = float(6.02214076e+23)
    mu = get_mu_col(df["Energy"], mu_df)
    print(mu.value_counts(dropna=False))
    # rho, d, mass, molar_mass, t_irr = float(rho), float(d), float(mass), float(molar_mass), float(t_irr)
    rho = float(rho)
    d = float(d)
    mu = mu.astype(np.float64)
    # print((mu*rho*d).max())
    print(f"mu is {mu}")
    print(f"rho is {rho}")
    print(f"d is {d}")
    k = (mu * rho * d) / (1 - np.exp(- mu * rho * d))
    lam = np.log(2) / df["Half-life [s]"]
    mass = float(mass)
    molar_mass = float(molar_mass)
    N = mass*AVOGADRO/molar_mass
    t_irr = float(t_irr)
    irr_start = pd.to_datetime(irr_start_str, dayfirst=True)
    print(f" irr start is {irr_start}")
    irr_end = irr_start + pd.to_timedelta(t_irr, "s")
    print(f"irr end is {irr_end}")
    print(f"acq. started on {df['Acquisition Started'][0]}")
    delta_t = pd.to_datetime(df["Acquisition Started"][0]) - irr_end
    delta_t = delta_t.total_seconds()
    delta_t = float(delta_t)
    print("foo pyco")
    print(f"lam is {lam}, type {type(lam)}")
    print(f"t_irr is {t_irr}, type {type(t_irr)}")
    print(f"delta_t is {delta_t}, type {type(delta_t)}")
    print(f"real time is {df['Real Time'][0]}, type {type(df['Real Time'][0])}")


    df["RR"] = (
                ((df['Real Time'][0] / df["Live Time"][0])
                 * k * lam * df["Area"])
                / (
                   N * (1-np.exp(-lam*t_irr)) * np.exp(-lam * delta_t)
                   * (1-np.exp(-lam * df["Real Time"][0]))
                   * df["eps"] * df["Ig"]))
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
