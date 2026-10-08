"""Euclid reference cosmology shared by Hi-COLA helper scripts.

The reference cosmology is matched to the Hi-COLA Lua parameter files:
h=0.67, Omega_cdm=0.270, Omega_b=0.049, Omega_Lambda=0.681,
n_s=0.96, A_s=2.10521e-9, and sigma8=0.816. The CMB temperature
and effective neutrino count are also matched so nbodykit/CLASS has a
nonzero radiation component.

Reference: Euclid Collaboration, "Euclid preparation: VII. Forecast
validation for Euclid cosmological probes", arXiv:1910.09273.
"""

from astropy import units as u
from astropy.cosmology import FlatLambdaCDM


SOURCE = "Euclid preparation: VII, arXiv:1910.09273"

h = 0.67
H0 = 100.0 * h * u.km / u.s / u.Mpc
Omega_cdm = 0.270
Ob0 = 0.049
Om0 = Omega_cdm + Ob0
Ode0 = 0.681
Tcmb0 = 2.7260 * u.K
Neff = 3.046
ns = 0.96
As = 2.10521e-9
sigma8 = 0.816
w0 = -1.0
wa = 0.0

PARAMETERS = {
    "H0": H0,
    "h": h,
    "Om0": Om0,
    "Omega_cdm": Omega_cdm,
    "Ob0": Ob0,
    "Ode0": Ode0,
    "Tcmb0": Tcmb0,
    "Neff": Neff,
    "ns": ns,
    "As": As,
    "sigma8": sigma8,
    "w0": w0,
    "wa": wa,
}

SERIALIZABLE_PARAMETERS = {
    "H0_km_s_Mpc": float(H0.to_value(u.km / u.s / u.Mpc)),
    "h": h,
    "Omega_m": Om0,
    "Omega_cdm": Omega_cdm,
    "Omega_b": Ob0,
    "Omega_Lambda": Ode0,
    "Tcmb0_K": float(Tcmb0.to_value(u.K)),
    "Neff": Neff,
    "n_s": ns,
    "A_s": As,
    "sigma8": sigma8,
    "w0": w0,
    "wa": wa,
}

FLATLAMBDACDM_KWARGS = {
    "H0": H0,
    "Om0": Om0,
    "Ob0": Ob0,
    "Tcmb0": Tcmb0,
    "Neff": Neff,
    "name": "Euclid reference cosmology",
}

EUCLID = FlatLambdaCDM(**FLATLAMBDACDM_KWARGS)


def get_cosmology() -> FlatLambdaCDM:
    """Return the Euclid reference cosmology for Astropy calculations."""
    return EUCLID
