#!/usr/bin/env python3

import argparse
import numpy as np
import MDAnalysis as mda
from tqdm import tqdm


def calculate_shortest_ion_projection(
    topology,
    trajectory,
    phosphate_selection,
    ion_selection,
    output,
):
    u = mda.Universe(topology, trajectory)

    r_arr = []
    phi_arr = []
    z_arr = []

    for ts in tqdm(u.trajectory, desc="Processing trajectory"):
        phosphates = u.select_atoms(phosphate_selection)
        ions = u.select_atoms(ion_selection)

        if len(phosphates) == 0:
            raise ValueError(f"No atoms found for phosphate selection: {phosphate_selection}")

        for ion in ions:
            distances = np.linalg.norm(phosphates.positions - ion.position, axis=1)
            nearest_idx = np.argmin(distances)

            nearest_p = phosphates[nearest_idx]

            dx, dy, dz = ion.position - nearest_p.position

            phi = np.arctan2(dy, dx)
            signed_r = np.sign(phi) * np.sqrt(dx**2 + dy**2)

            r_arr.append(signed_r)
            phi_arr.append(phi)
            z_arr.append(dz)

    output_data = np.column_stack((r_arr, phi_arr, z_arr))

    np.savetxt(
        output,
        output_data,
        fmt="%.4f",
        header="r_xy phi z",
        comments="# ",
    )

    print(f"Saved: {output}")


def main():
    parser = argparse.ArgumentParser(
        description="Project ions relative to their nearest selected atom, usually nearest Mg2+ to phosphate."
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
        "--phosphate-selection",
        default="name P",
        help="Reference atom selection. Default: 'name P'",
    )

    parser.add_argument(
        "--ion-selection",
        default="resname Mg",
        help="Ion atom selection. Default: 'resname Mg'",
    )

    parser.add_argument(
        "-o", "--output",
        default="2D-ion-projection-shortest-Mg.dat",
        help="Output filename. Default: 2D-ion-projection-shortest-Mg.dat",
    )

    args = parser.parse_args()

    calculate_shortest_ion_projection(
        topology=args.topology,
        trajectory=args.trajectory,
        phosphate_selection=args.phosphate_selection,
        ion_selection=args.ion_selection,
        output=args.output,
    )


if __name__ == "__main__":
    main()


# Example:
# python shortest_ion_projection.py \
#    -p rU30-Mg5-Na20-1.pdb \
#    -t principal-aligned.dcd \
#    --phosphate-selection "resname URA and name P" \
#    --ion-selection "resname Mg and around 10 resname URA" \
#    -o 2D-ion-projection-shortest-Mg.dat
