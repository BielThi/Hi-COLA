"""
GR vs nDGP: zoom on a filament, with the matter density in the background and
the haloes drawn as circles ("halo bubbles"), in the style of the f(R) movie.

Each panel shows
  * the projected matter density of a thin slab (smoothed, log scale),
  * a circle for every halo above M_MIN, with radius R200m (times CIRCLE_SCALE),
  * the log10 mass of the N_LABEL most massive haloes,
  * a header with the number of haloes above M_MIN and log10 of their total mass.

The runs must share initial conditions, so the same structures appear in both.
In nDGP you expect slightly more haloes and more mass around them.

Usage:
    python make_gr_ndgp_halos.py --gr-part GR_snap --gr-halo GR_halos \
                                 --mg-part NDGP_snap --mg-halo NDGP_halos --label "nDGP N1.2"
    python make_gr_ndgp_halos.py --demo              # synthetic data, to test the layout

The particle snapshots are optional. Without them the background is built from
the haloes themselves (mass-weighted and smoothed), which already looks like the movie.

Edit load_particles() and load_haloes() for your file formats.
Output: gr_ndgp_halos.pdf / .png / .svg (text kept as text).
"""
import argparse
import glob
import logging
import struct
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Circle
from scipy.ndimage import gaussian_filter
logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)

# ----------------------------------------------------------------------------
# Settings
# ----------------------------------------------------------------------------
BOX = 500.0              # box size in Mpc/h
OMEGA_M = 0.315192          # for R200m from the halo mass
CENTER = (0,0,0)            # (x, y, z) in Mpc/h; None = centre on the most massive GR halo
WIDTH = 500.0             # side of the shown region in Mpc/h (as in the movie)
DEPTH = 15.0             # slab thickness along the projection axis in Mpc/h
AXIS = 2                 # projection axis (0=x, 1=y, 2=z)
M_MIN = 10**12.5         # haloes drawn as circles above this mass (Msun/h)
N_LABEL = 2              # label the masses of the N most massive haloes
CIRCLE_SCALE = 1.0       # circle radius = CIRCLE_SCALE * R200m (increase if circles are too small)
NPIX = 500               # pixels per side of the background map
SMOOTH = 0.5             # smoothing of the background in Mpc/h
SCALE_BAR = 10.0         # scale bar in Mpc/h
FIG_W_IN = 6.7           # = \textwidth (17 cm)
DARK = True              # black background as in the movie; False for a white page
CMAP = "magma"
RES_LIM = 0.15           # --residual: colour scale of log10(rho_MG / rho_GR) is +-RES_LIM
RES_SMOOTH = 1.5         # --residual: smoothing in Mpc/h (small shifts of the same structures
                         # between the runs give red/blue pairs at the scale of SMOOTH)
# --residual colour map: black at 0, blue/cyan = less dense in MG, green/lime = denser
# (hues absent from magma, so not confused with the density panels)
RES_CMAP = LinearSegmentedColormap.from_list(
    "cyan_black_lime", ["#b5f0ff", "#2b8fd6", "#0b2a4a", "black", "#0f3d1c", "#3aa84a", "#d2f58a"])

SERIF = ["Times New Roman", "Liberation Serif", "DejaVu Serif"]
plt.rcParams.update({"font.family": "serif", "font.serif": SERIF,
                     "mathtext.fontset": "dejavuserif",
                     "pdf.fonttype": 42, "svg.fonttype": "none"})


# ----------------------------------------------------------------------------
# Reading: adapt to your files
# ----------------------------------------------------------------------------
def load_particles(path, center):
    """(N, 3) positions in Mpc/h of the particles in the slab around center
    (all 1024^3 positions would be 12 GB, so each file is read and cut separately).
      * Hi-COLA Gadget-2 snapshot (format 1, multi-file): give the prefix,
        e.g. .../snapshot_<name>_z0.500/gadget_z0.500 (reads gadget_z0.500.*)
      * numpy .npy with an (N, 3) array"""
    if path.endswith(".npy"):
        pos = np.load(path).astype(np.float32)
        return pos[in_slab(pos, center)[1]]
    files = sorted(glob.glob(path + ".*"), key=lambda f: int(f.rsplit(".", 1)[1]))
    if not files:
        raise FileNotFoundError(f"No Gadget files matched {path}.*")
    keep = []
    for f in files:
        with open(f, "rb") as fh:
            fh.seek(4)
            n = sum(struct.unpack("<6I", fh.read(24)))     # header: npart per type
            fh.seek(4 + 256 + 4 + 4)                       # skip header block, POS block marker
            pos = np.fromfile(fh, dtype="<f4", count=3 * n).reshape(n, 3)
        keep.append(pos[in_slab(pos, center)[1]])
    return np.concatenate(keep)


def load_haloes(path):
    """Return (pos (N,3) in Mpc/h, mass (N,) in Msun/h).
    Examples:
      * numpy .npz with arrays 'pos' and 'mass'
      * ROCKSTAR ascii:  d = np.loadtxt(path); pos = d[:, 8:11]; mass = d[:, 2]  (check the columns)
      * AHF / HODDIES halo files: read x, y, z and M200m
      * Hi-COLA halos200c_catalog.npy: structured array with x, y, z, m200c.
        Despite its name, m200c is the FoF mass (length * particle mass, b = 0.2),
        taken here as ~M200m to match R200m"""
    d = np.load(path)
    if path.endswith(".npy"):
        return np.c_[d["x"], d["y"], d["z"]].astype(np.float32), d["m200c"].astype(np.float64)
    return d["pos"].astype(np.float32), d["mass"].astype(np.float64)


# ----------------------------------------------------------------------------
def periodic_delta(x, c):
    return (x - c + BOX / 2) % BOX - BOX / 2


def in_slab(pos, center):
    """Projected coordinates (relative to the centre) of the points inside the slab."""
    d = periodic_delta(pos, np.asarray(center, np.float32))
    other = [i for i in range(3) if i != AXIS]
    m = np.abs(d[:, AXIS]) < DEPTH / 2
    for i in other:
        m &= np.abs(d[:, i]) < WIDTH / 2
    return d[:, other], m


def density_map(xy, weights=None, smooth=SMOOTH):
    h, _, _ = np.histogram2d(xy[:, 0], xy[:, 1], bins=NPIX, weights=weights,
                             range=[[-WIDTH / 2, WIDTH / 2]] * 2)
    h = gaussian_filter(h.T, smooth / (WIDTH / NPIX))
    return h / h[h > 0].mean()


def r200m(mass):
    rho_m = 2.775e11 * OMEGA_M            # mean matter density in Msun/h / (Mpc/h)^3
    return (3 * mass / (4 * np.pi * 200 * rho_m)) ** (1 / 3)


# ----------------------------------------------------------------------------
def demo_data():
    """Synthetic haloes along filaments (not physical). The nDGP copy has
    slightly more haloes and slightly higher masses, as expected."""
    rng = np.random.default_rng(3)
    nodes = rng.uniform(-20, 20, size=(6, 2))
    pts = []
    for a, b in zip(nodes[:-1], nodes[1:]):
        t = rng.uniform(0, 1, 70)[:, None]
        pts.append(a + t * (b - a) + rng.normal(0, 1.6, (70, 2)))
    xy = np.vstack(pts)
    pos = np.c_[xy + 250, rng.uniform(-6, 6, len(xy)) + 250].astype(np.float32)
    mass = 10 ** rng.normal(12.2, 0.45, len(xy))
    mass[np.argsort(mass)[-2:]] *= 8
    pos_mg = pos.copy()
    mass_mg = mass * 10 ** rng.normal(0.04, 0.03, len(mass))
    return (None, (pos, mass)), (None, (pos_mg, mass_mg))


def background(part, halo, center, smooth=SMOOTH):
    """log10 of the projected density (relative to its mean) in the slab."""
    hpos, hmass = halo
    hxy, hm = in_slab(hpos, center)
    if part is not None:
        pxy, pm = in_slab(part, center)
        bg = density_map(pxy[pm], smooth=smooth)
    else:                                        # background from the haloes themselves
        w = hmass[hm] ** (2 / 3)                 # compact blob per halo plus a diffuse filament
        bg = density_map(hxy[hm], weights=w, smooth=0.5) + \
             0.15 * density_map(hxy[hm], weights=w, smooth=3.0)
    return np.log10(np.maximum(bg, 1e-3))


def frame(ax, edge):
    """Scale bar (colour edge), limits, no ticks or spines."""
    x0 = -WIDTH / 2 + 0.06 * WIDTH
    ax.plot([x0, x0 + SCALE_BAR], [-WIDTH / 2 + 0.07 * WIDTH] * 2, color=edge, lw=1.4)
    ax.text(x0 + SCALE_BAR / 2, -WIDTH / 2 + 0.09 * WIDTH, rf"{SCALE_BAR:.0f} $h^{{-1}}$Mpc",
            ha="center", va="bottom", fontsize=6.5, color=edge)
    ax.set_xlim(-WIDTH / 2, WIDTH / 2); ax.set_ylim(-WIDTH / 2, WIDTH / 2)
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)


def panel(ax, bg, halo, center, title, vlim, sep=",  "):
    hpos, hmass = halo
    hxy, hm = in_slab(hpos, center)
    ext = [-WIDTH / 2, WIDTH / 2, -WIDTH / 2, WIDTH / 2]
    ax.imshow(bg, origin="lower", extent=ext, cmap=CMAP, vmin=vlim[0], vmax=vlim[1],
              interpolation="bilinear")

    sel = hm & (hmass > M_MIN)
    xy, m = hxy[sel], hmass[sel]
    edge = "white" if DARK else "black"
    for (x, y), mm in zip(xy, m):
        ax.add_patch(Circle((x, y), CIRCLE_SCALE * r200m(mm), fill=False,
                            ec=edge, lw=0.6, alpha=0.9))
    for i in np.argsort(m)[::-1][:N_LABEL]:
        r = CIRCLE_SCALE * r200m(m[i])
        ax.text(xy[i, 0] + 0.75 * r + 0.3, xy[i, 1] + 0.75 * r + 0.3,
                f"{np.log10(m[i]):.2f}", ha="left", va="bottom",
                fontsize=6.5, fontweight="bold", color=edge)
    head = (rf"{len(m)} haloes $>10^{{{np.log10(M_MIN):.1f}}}\,h^{{-1}}M_\odot${sep}"
            rf"$\log_{{10}}\Sigma M_{{\rm h}}={np.log10(m.sum()):.2f}$")
    ax.set_title(f"{title}\n{head}", fontsize=7.5, color=edge, pad=3, linespacing=1.4)
    frame(ax, edge)


def residual_panel(ax, diff, label):
    """log10(rho_MG / rho_GR) on a diverging scale of +-RES_LIM, with a colour bar."""
    edge = "white" if DARK else "black"
    ink = "white"                                # on the map, which is black around 0
    ext = [-WIDTH / 2, WIDTH / 2, -WIDTH / 2, WIDTH / 2]
    im = ax.imshow(diff, origin="lower", extent=ext, cmap=RES_CMAP, vmin=-RES_LIM, vmax=RES_LIM,
                   interpolation="bilinear")
    ax.set_title(f"{label} $-$ GR\n" + r"$\log_{10}(\rho_{\rm MG}/\rho_{\rm GR})$",
                 fontsize=7.5, color=edge, pad=3, linespacing=1.4)
    cax = ax.inset_axes([0.5, 0.07, 0.44, 0.03])
    cb = ax.figure.colorbar(im, cax=cax, orientation="horizontal", ticks=[-RES_LIM, 0, RES_LIM])
    cb.ax.tick_params(labelsize=6, colors=ink, length=2, pad=1)
    cb.outline.set_edgecolor(ink)
    frame(ax, ink)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gr-part"); ap.add_argument("--gr-halo")
    ap.add_argument("--mg-part"); ap.add_argument("--mg-halo")
    ap.add_argument("--label", default="nDGP N1.2")
    ap.add_argument("--z", default="0.5")
    ap.add_argument("--demo", action="store_true")
    ap.add_argument("--residual", action="store_true",
                    help="add a third panel with log10(rho_MG / rho_GR)")
    a = ap.parse_args()

    if a.demo:
        (gp, gh), (mp, mh) = demo_data()
    else:
        gh, mh = load_haloes(a.gr_halo), load_haloes(a.mg_halo)

    center = CENTER if CENTER is not None else gh[0][np.argmax(gh[1])]
    if a.demo:
        center = (250.0, 250.0, 250.0)
    elif a.gr_part or a.mg_part:                 # particles are cut to the slab while reading
        gp = load_particles(a.gr_part, center) if a.gr_part else None
        mp = load_particles(a.mg_part, center) if a.mg_part else None
    else:
        gp = mp = None
    vlim = (-1.0, 1.6)                       # same colour scale in both panels
    bgc = "black" if DARK else "white"

    n = 3 if a.residual else 2
    sep = "\n" if a.residual else ",  "      # narrower panels: header on two lines
    top = 0.76 if a.residual else 0.85
    fig, axes = plt.subplots(1, n, figsize=(FIG_W_IN, FIG_W_IN / n * 1.13 / top * 0.85),
                             facecolor=bgc)
    fig.subplots_adjust(left=0.01, right=0.99, top=top, bottom=0.01, wspace=0.03)
    bg_gr, bg_mg = background(gp, gh, center), background(mp, mh, center)
    panel(axes[0], bg_gr, gh, center, "GR", vlim, sep)
    panel(axes[1], bg_mg, mh, center, a.label, vlim, sep)
    if a.residual:
        residual_panel(axes[2], background(mp, mh, center, RES_SMOOTH)
                       - background(gp, gh, center, RES_SMOOTH), a.label)
    fig.suptitle(f"Same initial conditions, $z={a.z}$, field of view {WIDTH:.0f} $h^{{-1}}$Mpc",
                 fontsize=8, color="white" if DARK else "black", y=0.985)
    for ext_ in ("pdf", "png", "svg"):
        fig.savefig(f"gr_ndgp_halos.{ext_}", dpi=300, facecolor=bgc)
    print("wrote gr_ndgp_halos.pdf/.png/.svg")


if __name__ == "__main__":
    main()