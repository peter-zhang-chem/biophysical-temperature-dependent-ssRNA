#!/usr/bin/env python3

import argparse
import numpy as np
import MDAnalysis as md
from tqdm import tqdm


def calculate_ocf(topology, trajectory, selection="name P", start=1000, output="OCF.dat"):
    u = md.Universe(topology, trajectory)
    atoms = u.select_atoms(selection)

    n_atoms = len(atoms)
    n_bonds = n_atoms - 1

    if n_atoms < 2:
        raise ValueError(f"Selection '{selection}' contains fewer than 2 atoms.")

    ocf = np.zeros(n_bonds)
    distances = np.zeros(n_bonds)

    production = u.trajectory[start:]
    n_frames = len(production)

    if n_frames == 0:
        raise ValueError("No frames selected. Check the trajectory length and start frame.")

    for ts in tqdm(production, desc="Calculating OCF"):
        atoms = u.select_atoms(selection)

        vecs = []

        for i in range(n_bonds):
            vec = atoms[i].position - atoms[i + 1].position
            distance = np.linalg.norm(vec)

            if distance == 0:
                continue

            vec /= distance
            vecs.append(vec)
            distances[i] += distance

        for i in range(len(vecs)):
            for j in range(i + 1, len(vecs)):
                separation = j - i - 1
                ocf[separation] += np.dot(vecs[i], vecs[j])

    ocf /= n_frames * n_bonds
    distances /= n_frames
    avg_pp_distance = np.mean(distances)

    with open(output, "w") as f:
        f.write("# Orientation Correlation Function\n")
        f.write(f"# Topology: {topology}\n")
        f.write(f"# Trajectory: {trajectory}\n")
        f.write(f"# Atom selection: {selection}\n")
        f.write(f"# Start frame: {start}\n")
        f.write(f"# Average P-P distance: {avg_pp_distance:.6f}\n")
        f.write("# Separation_Index OCF\n")

        for i, value in enumerate(ocf):
            f.write(f"{i} {value:.8f}\n")

    print(f"Saved OCF data to {output}")
    print(f"Average P-P distance: {avg_pp_distance:.6f}")


def main():
    parser = argparse.ArgumentParser(
        description="Calculate the orientation correlation function from an MD trajectory."
    )

    parser.add_argument(
        "-p", "--topology",
        required=True,
        help="Topology file, e.g. system.pdb"
    )

    parser.add_argument(
        "-t", "--trajectory",
        required=True,
        help="Trajectory file, e.g. trajectory.dcd"
    )

    parser.add_argument(
        "-s", "--selection",
        default="name P",
        help="MDAnalysis atom selection. Default: 'name P'"
    )

    parser.add_argument(
        "--start",
        type=int,
        default=1000,
        help="Starting frame for analysis. Default: 1000"
    )

    parser.add_argument(
        "-o", "--output",
        default="OCF.dat",
        help="Output file name. Default: OCF.dat"
    )

    args = parser.parse_args()

    calculate_ocf(
        topology=args.topology,
        trajectory=args.trajectory,
        selection=args.selection,
        start=args.start,
        output=args.output,
    )


if __name__ == "__main__":
    main()


# Example:
# python calculate_ocf.py \
#    -p rA30-Mg10-Na300-1.pdb \
#    -t rA30-Mg10-Na300-1.dcd \
#    --start 1000 \
#    -o OCF.dat
