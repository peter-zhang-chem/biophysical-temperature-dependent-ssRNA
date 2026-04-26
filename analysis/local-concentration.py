#!/usr/bin/env python3

import argparse
import numpy as np
import MDAnalysis as mda
from tqdm import tqdm


def calculate_local_concentration(
    topology,
    trajectory,
    center_selection,
    ion_selection,
    output="local_concentration.dat",
    box_size=720.0,
    start=3000,
    end=-1,
    step=1,
    dr=0.1,
):
    u = mda.Universe(topology, trajectory)

    distances = []

    for ts in tqdm(u.trajectory[start:end:step], desc="Processing trajectory"):
        center_group = u.select_atoms(center_selection)
        ions = u.select_atoms(ion_selection)

        if len(center_group) == 0:
            raise ValueError(f"No atoms found for center selection: {center_selection}")

        if len(ions) == 0:
            continue

        center = center_group.center_of_geometry()
        d = np.linalg.norm(ions.positions - center, axis=1)
        distances.extend(d)

    frame_count = len(u.trajectory[start:end:step])

    if frame_count == 0:
        raise ValueError("No frames selected. Check --start, --end, and --step.")

    bin_edges = np.arange(0.1, box_size / 2, dr)
    hist, _ = np.histogram(distances, bins=bin_edges)

    r = bin_edges[:-1]
    shell_volumes = 4 * np.pi * r**2 * dr

    ion_count_avg = hist / frame_count
    local_conc = ion_count_avg / (shell_volumes * 6.022e-7)

    output_data = np.column_stack((r, local_conc))

    np.savetxt(
        output,
        output_data,
        fmt="%.6f",
        header="Distance Local_Concentration_mM",
        comments="",
    )

    print(f"Saved: {output}")


def main():
    parser = argparse.ArgumentParser(
        description="Calculate radial local ion concentration around the center of a selected molecule/group."
    )

    parser.add_argument(
        "-p", "--topology",
        required=True,
        help="Topology file, e.g. system.pdb",
    )

    parser.add_argument(
        "-t", "--trajectory",
        required=True,
        help="Trajectory file, e.g. trajectory.dcd",
    )

    parser.add_argument(
        "-c", "--center-selection",
        default="resname ADE",
        help="Selection used to define the radial origin. Default: 'resname ADE'",
    )

    parser.add_argument(
        "-i", "--ion-selection",
        default="resname Mg",
        help="Ion selection. Default: 'resname Mg'",
    )

    parser.add_argument(
        "-o", "--output",
        default="local_concentration.dat",
        help="Output filename. Default: local_concentration.dat",
    )

    parser.add_argument(
        "--box-size",
        type=float,
        required=True,
        help="Simulation box size in Angstrom.",
    )

    parser.add_argument(
        "--start",
        type=int,
        default=3000,
        help="Starting frame. Default: 3000",
    )

    parser.add_argument(
        "--end",
        type=int,
        default=-1,
        help="Ending frame. Default: -1",
    )

    parser.add_argument(
        "--step",
        type=int,
        default=1,
        help="Frame stride. Default: 1",
    )

    parser.add_argument(
        "--dr",
        type=float,
        default=0.1,
        help="Radial bin width in Angstrom. Default: 0.1",
    )

    args = parser.parse_args()

    calculate_local_concentration(
        topology=args.topology,
        trajectory=args.trajectory,
        center_selection=args.center_selection,
        ion_selection=args.ion_selection,
        output=args.output,
        box_size=args.box_size,
        start=args.start,
        end=args.end,
        step=args.step,
        dr=args.dr,
    )


if __name__ == "__main__":
    main()

# Example:
# python local_concentration.py \
#    -p rA30-Mg1-Na150-1.pdb \
#    -t md-align-wrap.dcd \
#    -c "resname ADE" \
#    -i "resname Mg" \
#    --box-size 720 \
#    --start 3000 \
#    -o local_concentration.dat
