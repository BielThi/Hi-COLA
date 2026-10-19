#!/usr/bin/env python3
"""Copy the measurements of one Hi-COLA run from scratch into $HICOLA_DATA_ROOT/hicola/.

For every snapshot_<simulation_name>_z<z>/ of the run it copies, when present:

    pofk_z<z>.txt                    -> HiCOLA_<simulation_name>_z<z>_PS.txt
    halo_catalog/halos200c_hmf.txt   -> HiCOLA_<simulation_name>_z<z>_HMF.txt   (+ .json)
    halo_catalog/hod_dev-v2/<name>_PS.txt   -> HiCOLA_<simulation_name>_z<z>_HOD-dev-v2_<name>_PS.txt   (+ .h5, .json)
    halo_catalog/hod_dev-v2/<name>_MPS.txt  -> HiCOLA_<simulation_name>_z<z>_HOD-dev-v2_<name>_MPS.txt  (+ .h5, .json)

Scratch is cleaned up automatically; these copies are the ones to keep. Called at the end of
jobs/backend.sh, jobs/halo_finder.sh and jobs/hod_mock.sh, and safe to rerun.

    python scripts/collect_products.py Abacus_Small nDGP_N1.2 s1234-fixed
"""

from __future__ import annotations

import argparse
import re
import shutil
import sys
from pathlib import Path

import hicola_paths as hp


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("target", help="e.g. COLA_NM4096, Abacus_Small")
    parser.add_argument("model")
    parser.add_argument("realization")
    args = parser.parse_args()
    run = (args.target, args.model, args.realization)

    run_dir = hp.run_dir(*run)
    prefix = f"snapshot_{hp.simulation_name(*run)}_z"
    snapshots = sorted(run_dir.glob(prefix + "*"))
    if not snapshots:
        print(f"No snapshot folders in {run_dir}", file=sys.stderr)
        return 1

    copied = 0
    for snapshot in snapshots:
        match = re.fullmatch(re.escape(prefix) + r"(\d+\.\d{3})", snapshot.name)
        if not match:
            continue
        redshift = float(match.group(1))
        sources = {
            "PS": [hp.pofk_file(*run, redshift)],
            "HMF": [hp.hmf_file(*run, redshift), hp.hmf_file(*run, redshift).with_suffix(".txt.json")],
        }
        for kind, files in sources.items():
            if not files[0].exists():
                continue
            destination = hp.hicola_file(*run, redshift, kind)
            destination.parent.mkdir(parents=True, exist_ok=True)
            for source in files:
                if source.exists():
                    target = destination.with_name(destination.name + source.name[len(files[0].name):])
                    shutil.copy2(source, target)
                    print(f"{source} -> {target}")
                    copied += 1
        hod_spectra = [
            (spectrum, statistic)
            for statistic in hp.HOD_STATISTICS
            for spectrum in sorted(hp.hod_dir(*run, redshift).glob(f"*_{statistic}.txt"))
        ]
        for spectrum, statistic in hod_spectra:
            name = spectrum.name[: -len(f"_{statistic}.txt")]
            destination = hp.hicola_hod_spectrum_file(*run, redshift, name, statistic)
            destination.parent.mkdir(parents=True, exist_ok=True)
            for source, target in [
                (spectrum, destination),
                (spectrum.with_suffix(".h5"), destination.with_suffix(".h5")),
                (Path(str(spectrum) + ".json"), Path(str(destination) + ".json")),
            ]:
                if source.exists():
                    shutil.copy2(source, target)
                    print(f"{source} -> {target}")
                    copied += 1

    print(f"Copied {copied} file(s) to {hp.hicola_dir(*run)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
