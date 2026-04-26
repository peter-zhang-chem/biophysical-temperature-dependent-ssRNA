#!/usr/bin/env python3

import argparse
import numpy as np
import MDAnalysis as mda
from tqdm import tqdm


RESNAME_MAP = {
    "A": "ADE",
    "U": "URA",
    "C": "CYT",
    "G": "GUA",
}


def count_inner_outer(
    topology,
    trajectory,
    residue_type,
    cutoff_inner=3.2,
    cutoff_outer=6.1,
    around_cutoff=20.0,
    start=0,
    end=-1,
    step=1,
    output="hist_data.dat",
):
    u = mda.Universe(topology, trajectory)

    if residue_type not in RESNAME_MAP:
        raise ValueError(f"Invalid residue type: {residue_type}. Choose from A, U, C, G.")

    resname = RESNAME_MAP[residue_type]

    inner_count = []
    outer_count = []

    for ts in tqdm(u.trajectory[start:end:step], desc="Processing trajectory"):
        Base = u.select_atoms(f"resname {resname} and name P")
        Mg = u.select_atoms(f"resname Mg and around {around_cutoff} resname {resname}")

        inner_tally = 0
        outer_tally = 0

        if len(Base) == 0 or len(Mg) == 0:
            inner_count.append(0)
            outer_count.append(0)
            continue

        for ion in Mg:
            distances = np.linalg.norm(Base.positions - ion.position, axis=1)
            shortest = np.min(distances)

            if shortest < cutoff_inner:
                inner_tally += 1
            elif cutoff_inner < shortest < cutoff_outer:
                outer_tally += 1

        inner_count.append(inner_tally)
        outer_count.append(outer_tally)

    # Save output
    with open(output, "w") as f:
        f.write("# inner_shell_bound   outer_shell_bound\n")
        for inner, outer in zip(inner_count, outer_count):
            f.write(f"{inner}   {outer}\n")

    print(f"Saved: {output}")


def main():
    parser = argparse.ArgumentParser(
        description="Count inner- and outer-shell Mg ions around RNA phosphates."
    )

    parser.add_argument(
        "-p", "--topology",
        required=True,
        help="Topology file (e.g., .pdb)",
    )

    parser.add_argument(
        "-t", "--trajectory",
        required=True,
        help="Trajectory file (e.g., .dcd)",
    )

    parser.add_argument(
        "-r", "--residue",
        required=True,
        choices=["A", "U", "C", "G"],
        help="RNA residue type (A, U, C, G)",
    )

    parser.add_argument(
        "-o", "--output",
        default="hist_data.dat",
        help="Output filename",
    )

    parser.add_argument(
        "--start",
        type=int,
        default=0,
        help="Starting frame",
    )

    parser.add_argument(
        "--end",
        type=int,
        default=-1,
        help="Ending frame",
    )

    parser.add_argument(
        "--step",
        type=int,
        default=1,
        help="Frame stride",
    )

    parser.add_argument(
        "--inner-cutoff",
        type=float,
        default=3.2,
        help="Inner-shell cutoff (Å)",
    )

    parser.add_argument(
        "--outer-cutoff",
        type=float,
        default=6.1,
        help="Outer-shell cutoff (Å)",
    )

    parser.add_argument(
        "--around",
        type=float,
        default=20.0,
        help="Distance cutoff for selecting nearby ions",
    )

    args = parser.parse_args()

    count_inner_outer(
        topology=args.topology,
        trajectory=args.trajectory,
        residue_type=args.residue,
        cutoff_inner=args.inner_cutoff,
        cutoff_outer=args.outer_cutoff,
        around_cutoff=args.around,
        start=args.start,
        end=args.end,
        step=args.step,
        output=args.output,
    )


if __name__ == "__main__":
    main()

# Example
# python count_inner_outer.py \
#     -p rC30-Mg5-Na20-1.pdb \
#    -t md-align-wrap.dcd \
#    -r C \
#    --start 3000 \
#    -o hist_data.dat
