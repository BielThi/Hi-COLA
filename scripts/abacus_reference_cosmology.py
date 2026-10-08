"""AbacusSummit c000 (Planck18) reference cosmology shared by Hi-COLA helper scripts.

The reference cosmology is matched to the Hi-COLA Lua parameter files for the Abacus-flavored
runs (nDGP_Abacus.lua, LCDM_Abacus.lua, nDGP1dot0_Abacus.lua): h=0.6736, Omega_cdm=0.265890,
Omega_b=0.049302, Omega_Lambda=0.684808, n_s=0.9649, A_s=2.083e-9, and sigma8_cb=0.811355.
Omega_cdm here already includes AbacusSummit c000's massive-neutrino density (physically
Omega_nu=0.001420), folded in rather than kept separate: with
ic_type_of_input="powerspectrum", a nonzero cosmology_OmegaMNu makes the COLA solver turn on
scaledependent growth but never build the transfer data it then needs, crashing at startup.
Folding the neutrino mass into OmegaCDM instead keeps Omega_m matched to AbacusSummit c000
(0.315192) and avoids that path -- see the comment in nDGP_Abacus.lua for the full
explanation. The CMB temperature and effective neutrino count are also matched so
nbodykit/CLASS has a nonzero radiation component.

Reference: AbacusSummit cosmology c000, the Planck2018
base_plikHM_TTTEEE_lowl_lowE_lensing mean.
https://abacussummit.readthedocs.io/en/latest/cosmologies.html
"""

from astropy import units as u
from astropy.cosmology import FlatLambdaCDM


h = 0.6736
H0 = 100.0 * h * u.km / u.s / u.Mpc
Omega_cdm = 0.265890
Ob0 = 0.049302
Om0 = Omega_cdm + Ob0
Ode0 = 0.684808
Tcmb0 = 2.7260 * u.K
Neff = 3.046
ns = 0.9649
As = 2.083e-9
sigma8 = 0.811355
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
    "name": "AbacusSummit c000 (Planck18) reference cosmology",
}

ABACUS = FlatLambdaCDM(**FLATLAMBDACDM_KWARGS)


def get_cosmology() -> FlatLambdaCDM:
    """Return the AbacusSummit c000 (Planck18) reference cosmology for Astropy calculations."""
    return ABACUS
