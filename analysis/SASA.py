#!/usr/bin/env python3

import argparse
import re
import numpy as np
import MDAnalysis as mda
import freesasa
from tqdm import tqdm


class CGClassifier(freesasa.Classifier):
    purePython = True

    def classify(self, residueName, atomName):
        a = atomName.strip()

        if re.match(r"^P3\b", a) or re.match(r"^P\b", a):
            return "Phosphate"
        if re.match(r"^S\b", a):
            return "Sugar"
        if re.match(r"^(A|U|C|G)\b", a):
            return "Base"

        return "Other"

    def radius(self, residueName, atomName):
        a = atomName.strip()

        if re.match(r"^P3\b", a) or re.match(r"^P\b", a):
            return 1.89
        if re.match(r"^S\b", a):
            return 2.61
        if re.match(r"^A\b", a):
            return 2.52
        if re.match(r"^U\b", a):
            return 2.43
        if re.match(r"^C\b", a):
            return 2.41
        if re.match(r"^G\b", a):
            return 2.73

        return 0.0


def sasa_partitioned(u, classifier, atom_selection, probe_radius=2.0, n_slices=100):
    params = freesasa.Parameters(
        {
            "algorithm": freesasa.LeeRichards,
            "n-slices": int(n_slices),
            "probe-radius": float(probe_radius),
        }
    )

    atoms = u.select_atoms(atom_selection)

    if len(atoms) == 0:
        raise ValueError(f"No atoms found for selection: {atom_selection}")

    coords = atoms.positions.astype(float)
    resnames = atoms.resnames
    atom_names = atoms.names

    radii = np.array(
        [classifier.radius(resnames[i], atom_names[i]) for i in range(atoms.n_atoms)],
        dtype=float,
    )

    result = freesasa.calcCoord(coords.reshape(-1), radii, params)

    areas = {
        "Base": 0.0,
        "Sugar": 0.0,
        "Phosphate": 0.0,
        "Total": result.totalArea(),
    }

    for i in range(atoms.n_atoms):
        cls = classifier.classify(resnames[i], atom_names[i])
        if cls in areas:
            areas[cls] += result.atomArea(i)

    return areas


def calculate_sasa(
    topology,
    trajectory,
    atom_selection,
    output,
    start=0,
    end=None,
    step=1,
    probe_radius=2.0,
    n_slices=100,
):
    u = mda.Universe(topology, trajectory)
    classifier = CGClassifier()

    with open(output, "w") as f:
        f.write("# Frame Base Sugar Phosphate Total\n")

        for ts in tqdm(u.trajectory[start:end:step], desc="Calculating SASA"):
            areas = sasa_partitioned(
                u,
                classifier,
                atom_selection=atom_selection,
                probe_radius=probe_radius,
                n_slices=n_slices,
            )

            f.write(
                f"{ts.frame:8d} "
                f"{areas['Base']:12.6f} "
                f"{areas['Sugar']:12.6f} "
                f"{areas['Phosphate']:12.6f} "
                f"{areas['Total']:12.6f}\n"
            )

    print(f"Saved SASA data to {output}")


def main():
    parser = argparse.ArgumentParser(
        description="Calculate partitioned SASA for coarse-grained RNA beads using FreeSASA."
    )

    parser.add_argument(
        "-p", "--topology",
        required=True,
        help="Topology file, e.g. system.pdb",
    )

    parser.add_argument(
        "-t", "--trajectory",
        required=True,
        help="Trajectory file, e.g. traj.dcd",
    )

    parser.add_argument(
        "-s", "--selection",
        default="name A U C G S P P3",
        help="Atom selection for SASA calculation. Default: 'name A U C G S P P3'",
    )

    parser.add_argument(
        "-o", "--output",
        default="sasa_partitioned_cg.dat",
        help="Output filename. Default: sasa_partitioned_cg.dat",
    )

    parser.add_argument(
        "--start",
        type=int,
        default=0,
        help="Starting frame. Default: 0",
    )

    parser.add_argument(
        "--end",
        type=int,
        default=None,
        help="Ending frame. Default: use all frames after start",
    )

    parser.add_argument(
        "--step",
        type=int,
        default=10,
        help="Frame stride. Default: 10",
    )

    parser.add_argument(
        "--probe-radius",
        type=float,
        default=2.0,
        help="Probe radius in Å. Default: 2.0",
    )

    parser.add_argument(
        "--n-slices",
        type=int,
        default=100,
        help="Number of Lee-Richards slices. Default: 100",
    )

    args = parser.parse_args()

    calculate_sasa(
        topology=args.topology,
        trajectory=args.trajectory,
        atom_selection=args.selection,
        output=args.output,
        start=args.start,
        end=args.end,
        step=args.step,
        probe_radius=args.probe_radius,
        n_slices=args.n_slices,
    )


if __name__ == "__main__":
    main()

# Example
# python calculate_cg_sasa.py \
#    -p rU30-Mg5-Na20-1.pdb \
#    -t md2.dcd \
#    -s "name A U S P P3" \
#    --step 10 \
#    -o sasa_partitioned_cg.dat
