"""
GR vs nDGP: zoom on a filament, with the matter density in the background and
the haloes drawn as circles ("halo bubbles"), in the style of the f(R) movie.

Panels
  (a) GR and (b) nDGP: projected matter density of a thin slab (log scale), a circle
      for every halo above M_MIN (radius R200m x CIRCLE_SCALE), the log10 mass of the
      N_LABEL most massive haloes, and a header with the number and total mass of haloes.
      In (b) each circle is coloured by the mass change of that halo relative to its
      GR partner (red = heavier in nDGP, blue = lighter).
  (c) ln(rho_nDGP / rho_GR): what gravity changes once the shared structure is removed.

The script also prints box-wide numbers for the text: the median mass change of
matched haloes above M_MIN and the fraction that are heavier in nDGP.

The runs must share initial conditions, so the same structures appear in both.
In nDGP you expect slightly more haloes and more mass around them.

Usage:
    python make_gr_ndgp_halos2.py --gr-part GR_snap --gr-halo GR_halos \
                                  --mg-part NDGP_snap --mg-halo NDGP_halos --label "nDGP N1.2"
    python make_gr_ndgp_halos2.py ... --center-rank 20   # centre on the 20th most massive GR halo
    python make_gr_ndgp_halos2.py ... --center 120 340 55 # or give the centre in Mpc/h
    python make_gr_ndgp_halos2.py --demo              # synthetic data, to test the layout

Tip: do not centre on the most massive cluster. Vainshtein screening is strongest there;
the fifth force acts in filaments and voids. A halo of ~1e14 Msun/h on a filament is better.

The particle snapshots are optional. Without them the background is built from
the haloes themselves (mass-weighted and smoothed), which already looks like the movie.

Edit load_particles() and load_haloes() for your file formats.
Output: gr_ndgp_halos2.pdf / .png / .svg (text kept as text).
"""
import argparse
import glob
import logging
import struct
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from scipy.ndimage import gaussian_filter
from scipy.spatial import cKDTree
from matplotlib.colors import LinearSegmentedColormap, Normalize
import matplotlib.cm as cm
logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)

# ----------------------------------------------------------------------------
# Settings
# ----------------------------------------------------------------------------
BOX = 500.0              # box size in Mpc/h
OMEGA_M = 0.315192       # for R200m from the halo mass
CENTER = None            # (x, y, z) in Mpc/h; None = centre on the most massive GR halo
WIDTH = 46.0             # side of the shown region in Mpc/h (as in the movie)
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
RATIO_SMOOTH = 1.5       # extra smoothing for the ratio panel in Mpc/h (the ratio is noisy on small scales)
RATIO_LIM = 0.35         # colour range of ln(rho_nDGP/rho_GR), i.e. about +-40 % (= +-0.15 dex)
DLOGM_LIM = 0.08         # colour range of log10(M_nDGP/M_GR) for the circles (0.08 dex = 20 %)
MATCH_FRAC = 1.0         # a GR and an nDGP halo are partners if closer than MATCH_FRAC * R200m(GR)
                         # (Hi-COLA FoF centres move by ~0.7 R200m between the runs: 0.5 matches 31 %)
MATCH_MAX_DEX = 0.3      # ... and their masses differ by less than this (rejects wrong matches)
CMAP_RATIO = "RdBu_r"    # halo circles in (b): white at 0, visible on the dark density map
# ratio panel (c): black at 0, blue/cyan = emptier in nDGP, green/lime = denser
# (hues absent from magma, so not confused with the density panels)
CMAP_DIFF = LinearSegmentedColormap.from_list(
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


def match_haloes(gr, mg, m_min=M_MIN):
    """Pair every GR halo above m_min with the nearest nDGP halo (periodic box).
    Returns indices (i_gr, i_mg) of accepted pairs. Works because the runs share
    initial conditions, so the same halo sits at almost the same place."""
    (pg, mgr), (pm, mmg) = gr, mg
    sel_g = np.where(mgr > m_min)[0]
    sel_m = np.where(mmg > m_min / 3)[0]          # allow partners just below the threshold
    tree = cKDTree(np.mod(pm[sel_m], BOX), boxsize=BOX)
    d, j = tree.query(np.mod(pg[sel_g], BOX))
    i_mg = sel_m[j]
    ok = (d < MATCH_FRAC * r200m(mgr[sel_g])) & \
         (np.abs(np.log10(mmg[i_mg] / mgr[sel_g])) < MATCH_MAX_DEX)
    return sel_g[ok], i_mg[ok]


def report_matches(gr, mg):
    ig, im = match_haloes(gr, mg)
    dlog = np.log10(mg[1][im] / gr[1][ig])
    n_g = np.sum(gr[1] > M_MIN); n_m = np.sum(mg[1] > M_MIN)
    print(f"Box-wide: {n_g} GR and {n_m} nDGP haloes above {M_MIN:.2e} Msun/h "
          f"({100 * (n_m / n_g - 1):+.1f} %)")
    print(f"  matched pairs: {len(ig)} ({100 * len(ig) / n_g:.0f} % of GR haloes)")
    print(f"  median mass change: {100 * (10**np.median(dlog) - 1):+.1f} %, "
          f"heavier in nDGP: {100 * np.mean(dlog > 0):.0f} %")
    for lo, hi in [(12.5, 13), (13, 13.5), (13.5, 14), (14, 16)]:
        s = (np.log10(gr[1][ig]) >= lo) & (np.log10(gr[1][ig]) < hi)
        if s.sum() > 5:
            print(f"    log M in [{lo}, {hi}): N={s.sum():6d}  median dM/M = "
                  f"{100 * (10**np.median(dlog[s]) - 1):+.1f} %")
    return ig, im


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


def background(part, halo, center):
    hpos, hmass = halo
    hxy, hm = in_slab(hpos, center)
    if part is not None:
        pxy, pm = in_slab(part, center)
        return density_map(pxy[pm])
    w = hmass[hm] ** (2 / 3)                     # compact blob per halo plus a diffuse filament
    return density_map(hxy[hm], weights=w, smooth=0.5) + \
        0.15 * density_map(hxy[hm], weights=w, smooth=3.0)


def panel(ax, part, halo, center, title, vlim, dlogm=None):
    """dlogm: optional array (same length as halo) with log10(M_nDGP/M_GR) of the
    matched partner, NaN where unmatched; used to colour the circles."""
    hpos, hmass = halo
    hxy, hm = in_slab(hpos, center)
    if part is not None:
        pxy, pm = in_slab(part, center)
        bg = density_map(pxy[pm])
    else:                                        # background from the haloes themselves
        w = hmass[hm] ** (2 / 3)                 # compact blob per halo plus a diffuse filament
        bg = density_map(hxy[hm], weights=w, smooth=0.5) + \
             0.15 * density_map(hxy[hm], weights=w, smooth=3.0)
    bg = np.log10(np.maximum(bg, 1e-3))
    ext = [-WIDTH / 2, WIDTH / 2, -WIDTH / 2, WIDTH / 2]
    ax.imshow(bg, origin="lower", extent=ext, cmap=CMAP, vmin=vlim[0], vmax=vlim[1],
              interpolation="bilinear")

    sel = hm & (hmass > M_MIN)
    xy, m = hxy[sel], hmass[sel]
    edge = "white" if DARK else "black"
    dl = dlogm[sel] if dlogm is not None else np.full(len(m), np.nan)
    cmap_c = plt.get_cmap(CMAP_RATIO); norm_c = Normalize(-DLOGM_LIM, DLOGM_LIM)
    for (x, y), mm, d in zip(xy, m, dl):
        c = edge if np.isnan(d) else cmap_c(norm_c(d))
        ax.add_patch(Circle((x, y), CIRCLE_SCALE * r200m(mm), fill=False,
                            ec=c, lw=0.6 if np.isnan(d) else 1.1, alpha=0.95))
    for i in np.argsort(m)[::-1][:N_LABEL]:
        r = CIRCLE_SCALE * r200m(m[i])
        ax.text(xy[i, 0] + 0.75 * r + 0.3, xy[i, 1] + 0.75 * r + 0.3,
                f"{np.log10(m[i]):.2f}", ha="left", va="bottom",
                fontsize=6.5, fontweight="bold", color=edge)
    head = rf"{len(m)} haloes, $\log_{{10}}\Sigma M_{{\rm h}}={np.log10(m.sum()):.2f}$"
    ax.set_title(f"{title}\n{head}", fontsize=7.5, color=edge, pad=3, linespacing=1.4)

    x0 = -WIDTH / 2 + 0.06 * WIDTH
    ax.plot([x0, x0 + SCALE_BAR], [-WIDTH / 2 + 0.07 * WIDTH] * 2, color=edge, lw=1.4)
    ax.text(x0 + SCALE_BAR / 2, -WIDTH / 2 + 0.09 * WIDTH, rf"{SCALE_BAR:.0f} $h^{{-1}}$Mpc",
            ha="center", va="bottom", fontsize=6.5, color=edge)
    ax.set_xlim(ext[:2]); ax.set_ylim(ext[2:])
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gr-part"); ap.add_argument("--gr-halo")
    ap.add_argument("--mg-part"); ap.add_argument("--mg-halo")
    ap.add_argument("--label", default="nDGP N1.2")
    ap.add_argument("--z", default="0.5")
    ap.add_argument("--demo", action="store_true")
    ap.add_argument("--center", type=float, nargs=3, help="centre (x y z) in Mpc/h")
    ap.add_argument("--center-rank", type=int, default=1,
                    help="centre on the n-th most massive GR halo (default 1)")
    a = ap.parse_args()

    if a.demo:
        (gp, gh), (mp, mh) = demo_data()
    else:
        gh, mh = load_haloes(a.gr_halo), load_haloes(a.mg_halo)

    if a.center:
        center = tuple(a.center)
    elif CENTER is not None:
        center = CENTER
    else:
        center = gh[0][np.argsort(gh[1])[::-1][a.center_rank - 1]]
    if a.demo:
        center = (250.0, 250.0, 250.0)
    else:                                    # particles are cut to the slab while reading
        gp = load_particles(a.gr_part, center) if a.gr_part else None
        mp = load_particles(a.mg_part, center) if a.mg_part else None
    print("centre:", np.round(center, 2))

    # halo matching and the mass change of every nDGP halo relative to its GR partner
    ig, im = report_matches(gh, mh)
    dlog_mg = np.full(len(mh[1]), np.nan)
    dlog_mg[im] = np.log10(mh[1][im] / gh[1][ig])

    bg_gr, bg_mg = background(gp, gh, center), background(mp, mh, center)
    k = RATIO_SMOOTH / (WIDTH / NPIX)
    eps = 1e-2                               # avoids 0/0 in empty pixels
    ratio = np.log((gaussian_filter(bg_mg, k) + eps) / (gaussian_filter(bg_gr, k) + eps))
    if gp is None:
        print("note: without particles the ratio panel is built from haloes only")

    vlim = (-1.0, 1.6)                       # same colour scale in both density panels
    bgc = "black" if DARK else "white"
    fg = "white" if DARK else "black"
    fig, axes = plt.subplots(1, 3, figsize=(FIG_W_IN, FIG_W_IN / 3 * 1.42), facecolor=bgc)
    fig.subplots_adjust(left=0.01, right=0.99, top=0.80, bottom=0.17, wspace=0.03)
    panel(axes[0], gp, gh, center, "GR", vlim)
    panel(axes[1], mp, mh, center, a.label, vlim, dlogm=dlog_mg)

    ext = [-WIDTH / 2, WIDTH / 2, -WIDTH / 2, WIDTH / 2]
    ax = axes[2]
    imr = ax.imshow(ratio, origin="lower", extent=ext, cmap=CMAP_DIFF,
                    vmin=-RATIO_LIM, vmax=RATIO_LIM, interpolation="bilinear")
    ax.set_title("Difference\n" + r"$\ln(\rho_{\rm nDGP}/\rho_{\rm GR})$", fontsize=7.5,
                 color=fg, pad=3, linespacing=1.4)
    x0 = -WIDTH / 2 + 0.06 * WIDTH                     # scale bar, white on the dark map
    ax.plot([x0, x0 + SCALE_BAR], [-WIDTH / 2 + 0.07 * WIDTH] * 2, color="white", lw=1.4)
    ax.text(x0 + SCALE_BAR / 2, -WIDTH / 2 + 0.09 * WIDTH, rf"{SCALE_BAR:.0f} $h^{{-1}}$Mpc",
            ha="center", va="bottom", fontsize=6.5, color="white")
    ax.set_xlim(ext[:2]); ax.set_ylim(ext[2:])
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)

    cax1 = fig.add_axes([0.355, 0.10, 0.27, 0.025])
    cb1 = fig.colorbar(cm.ScalarMappable(Normalize(-DLOGM_LIM, DLOGM_LIM), CMAP_RATIO),
                       cax=cax1, orientation="horizontal",
                       ticks=[-DLOGM_LIM, -DLOGM_LIM / 2, 0, DLOGM_LIM / 2, DLOGM_LIM])
    cb1.set_label(r"halo mass change $\log_{10}(M_{\rm nDGP}/M_{\rm GR})$", fontsize=6.5,
                  color=fg, labelpad=1)
    cax2 = fig.add_axes([0.69, 0.10, 0.27, 0.025])
    cb2 = fig.colorbar(imr, cax=cax2, orientation="horizontal")
    cb2.set_label("denser (green) or emptier (blue) in nDGP", fontsize=6.5, color=fg, labelpad=1)
    for c in (cb1, cb2):
        c.ax.tick_params(labelsize=6, length=2, colors=fg)
        c.outline.set_edgecolor(fg)
    fig.suptitle(f"Same initial conditions, $z={a.z}$, field of view {WIDTH:.0f} $h^{{-1}}$Mpc, "
                 rf"haloes above $10^{{{np.log10(M_MIN):.1f}}}\,h^{{-1}}M_\odot$",
                 fontsize=8, color=fg, y=0.985)
    for ext_ in ("pdf", "png", "svg"):
        fig.savefig(f"gr_ndgp_halos2.{ext_}", dpi=300, facecolor=bgc)
    print("wrote gr_ndgp_halos2.pdf/.png/.svg")


if __name__ == "__main__":
    main()