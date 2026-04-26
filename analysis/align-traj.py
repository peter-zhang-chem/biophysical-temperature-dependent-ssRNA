#!/usr/bin/env python3

import argparse
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis import transformations as trans


def align_and_wrap_trajectory(
    topology,
    trajectory,
    selection,
    output="md-align-wrap.dcd",
    box_size=500.0,
    center="geometry",
):
    u = mda.Universe(topology, trajectory)

    group = u.select_atoms(selection)

    if len(group) == 0:
        raise ValueError(f"No atoms found for selection: {selection}")

    dim = np.array([box_size, box_size, box_size, 90, 90, 90])

    workflow = [
        trans.boxdimensions.set_dimensions(dim),
        trans.unwrap(u.atoms),
        trans.center_in_box(group, center=center),
        trans.wrap(u.atoms),
    ]

    u.trajectory.add_transformations(*workflow)

    align.AlignTraj(
        u,
        u,
        select=selection,
        filename=output,
        match_atoms=True,
    ).run()

    print(f"Saved aligned trajectory: {output}")


def main():
    parser = argparse.ArgumentParser(
        description="Unwrap, center, wrap, and align an MD trajectory using MDAnalysis."
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
        "-s", "--selection",
        required=True,
        help="Atom selection used for centering and alignment, e.g. 'resname GUA'",
    )

    parser.add_argument(
        "-o", "--output",
        default="md-align-wrap.dcd",
        help="Output aligned trajectory. Default: md-align-wrap.dcd",
    )

    parser.add_argument(
        "--box-size",
        type=float,
        required=True,
        help="Cubic box size in Å.",
    )

    parser.add_argument(
        "--center",
        choices=["geometry", "mass"],
        default="geometry",
        help="Centering method. Default: geometry",
    )

    args = parser.parse_args()

    align_and_wrap_trajectory(
        topology=args.topology,
        trajectory=args.trajectory,
        selection=args.selection,
        output=args.output,
        box_size=args.box_size,
        center=args.center,
    )


if __name__ == "__main__":
    main()


# Example
# python align_wrap_trajectory.py \
#    -p A-Mg5-Na20-1.pdb \
#    -t A-Mg5-Na20-1.dcd \
#    -s "resname ADE" \
#    --box-size 421 \
#    -o md-align-wrap.dcd
