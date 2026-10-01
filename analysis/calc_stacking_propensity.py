import argparse
import numpy as np
import MDAnalysis as mda
from MDAnalysis.lib.distances import calc_dihedrals
import MDAnalysis.transformations as trans


# dihedral angle definitions
PHI10_EQUIL = -2.58684
PHI20_EQUIL = 3.07135
STACK_TOLERANCE = 0.25
MIN_STACK_LENGTH = 3


def find_consecutive_stack_lengths(stack_states, min_length=3):
    """
    Find lengths of consecutive stacked ('S') segments.

    Parameters
    ----------
    stack_states : list[str]
        Sequence containing "S" and "N".
    min_length : int
        Minimum consecutive stack length to retain.

    Returns
    -------
    list[int]
        Lengths of consecutive stacked segments.
    """
    lengths = []
    count = 0

    for state in stack_states:
        if state == "S":
            count += 1
        else:
            if count >= min_length:
                lengths.append(count)
            count = 0

    # Handle a stacked segment that reaches the end
    if count >= min_length:
        lengths.append(count)

    return lengths


def is_stacked(phi10, phi20):
    """
    Determine whether a nucleotide step satisfies
    the stacking-angle criterion.
    """

    phi10_condition = (
        np.abs(phi10 - PHI10_EQUIL)
        / np.abs(PHI10_EQUIL)
        < STACK_TOLERANCE
    )

    phi20_condition = (
        np.abs(phi20 - PHI20_EQUIL)
        / np.abs(PHI20_EQUIL)
        < STACK_TOLERANCE
    )

    return phi10_condition and phi20_condition


def calculate_stack_states(P, S, B):
    """
    Calculate stacking states for one trajectory frame.

    Returns
    -------
    list[str]
        "S" for stacked and "N" for non-stacked.
    """

    phi10 = []
    phi20 = []

    # Internal nucleotide steps
    for i in range(len(P) - 2):

        phi10.append(
            calc_dihedrals(
                P[i].position,
                S[i].position,
                P[i + 1].position,
                S[i + 1].position,
            )
        )

        phi20.append(
            calc_dihedrals(
                S[i].position,
                P[i + 1].position,
                S[i + 1].position,
                P[i + 2].position,
            )
        )

    # Last nucleotide step
    phi10.append(
        calc_dihedrals(
            P[-2].position,
            S[-2].position,
            P[-1].position,
            S[-1].position,
        )
    )

    phi20.append(
        calc_dihedrals(
            S[-2].position,
            P[-1].position,
            S[-1].position,
            B[-1].position,
        )
    )

    stack_states = []

    for angle10, angle20 in zip(phi10, phi20):

        if is_stacked(angle10, angle20):
            stack_states.append("S")
        else:
            stack_states.append("N")

    return stack_states


# Wrapping and translate RNA to the center of the box
def add_wrapping_transformations(
    u,
    rna_selection,
    box_size,
):
    """
    Apply PBC transformations:

    1. Set cubic box dimensions
    2. Unwrap atoms
    3. Center RNA in the box
    4. Wrap all atoms back into the box

    Parameters
    ----------
    u : MDAnalysis.Universe
        MDAnalysis Universe.
    rna_selection : str
        Atom selection used to define the RNA for centering.
        Example: "resname ADE"
    box_size : float
        Cubic box length in Angstrom.
    """

    RNA = u.select_atoms(rna_selection)

    if len(RNA) == 0:
        raise ValueError(
            f"No atoms found for RNA selection: {rna_selection}"
        )

    dimensions = [
        box_size,
        box_size,
        box_size,
        90,
        90,
        90,
    ]

    workflow = [
        trans.boxdimensions.set_dimensions(dimensions),
        trans.unwrap(u.atoms),
        trans.center_in_box(
            RNA,
            center="mass",
        ),
        trans.wrap(u.atoms),
    ]

    u.trajectory.add_transformations(*workflow)

    print("Applied trajectory wrapping/centering transformations.")
    print(f"RNA selection: {rna_selection}")
    print(f"Box size:      {box_size:.2f} Å")


def analyze_stacking(
    topology,
    trajectory,
    start_frame=1000,
    stop_frame=None,
    step=1,
    nucleotide_name="A",
    output="stack_propensity.dat",
    apply_wrap=False,
    box_size=None,
    rna_selection=None,
):
    """
    Calculate RNA stacking propensity.

    Parameters
    ----------
    topology : str
        Topology/PDB file.
    trajectory : str
        Trajectory file.
    start_frame : int
        First trajectory frame.
    stop_frame : int or None
        Final trajectory frame.
    step : int
        Frame stride.
    nucleotide_name : str
        Base bead atom name, e.g. A, U, C, G.
    output : str
        Output filename.
    apply_wrap : bool
        Whether to apply PBC wrapping/centering.
    box_size : float or None
        Cubic box size in Angstrom.
    rna_selection : str or None
        MDAnalysis selection used to center RNA.
    """

    print(f"MDAnalysis version: {mda.__version__}")

    u = mda.Universe(
        topology,
        trajectory,
    )

    # Trajectory wrapping
    if apply_wrap:

        if box_size is None:
            raise ValueError(
                "--box-size must be provided when --wrap is used."
            )

        if rna_selection is None:
            raise ValueError(
                "--rna-selection must be provided when --wrap is used."
            )

        add_wrapping_transformations(
            u=u,
            rna_selection=rna_selection,
            box_size=box_size,
        )

    # Backbond selection
    P = u.select_atoms("name P")
    S = u.select_atoms("name S")
    B = u.select_atoms(f"name {nucleotide_name}")

    n_residues = len(P)

    if n_residues == 0:
        raise ValueError("No atoms with name P were found.")

    if len(S) == 0:
        raise ValueError("No atoms with name S were found.")

    if len(B) == 0:
        raise ValueError(
            f"No atoms with name {nucleotide_name} were found."
        )

    print()
    print("Analysis setup")
    print("--------------")
    print(f"Topology:        {topology}")
    print(f"Trajectory:      {trajectory}")
    print(f"RNA length:      {n_residues}")
    print(f"Total frames:    {len(u.trajectory)}")
    print(f"Start frame:     {start_frame}")
    print(f"Stop frame:      {stop_frame}")
    print(f"Frame step:      {step}")
    print(f"Wrapping:        {apply_wrap}")
    print()

    # Analysis begin
    results = []

    for ts in u.trajectory[start_frame:stop_frame:step]:

        stack_states = calculate_stack_states(
            P,
            S,
            B,
        )

        stack_lengths = find_consecutive_stack_lengths(
            stack_states,
            min_length=MIN_STACK_LENGTH,
        )

        # stacking propensity definition
        stacking_propensity = (
            np.sum(stack_lengths) / n_residues
        )

        if stack_lengths:
            n_max = np.max(stack_lengths)
            n_min = np.min(stack_lengths)
        else:
            n_max = 0
            n_min = 0

        results.append(
            [
                ts.frame,
                stacking_propensity,
                n_max,
                n_min,
            ]
        )

    results = np.asarray(results)

    np.savetxt(
        output,
        results,
        fmt=["%d", "%.4f", "%d", "%d"],
        header=(
            "frame "
            "stacking_propensity "
            "N_max "
            "N_min"
        ),
    )

    print(f"Saved results to: {output}")


def main():

    parser = argparse.ArgumentParser(
        description=(
            "Calculate RNA stacking propensity from an MD trajectory."
        )
    )

    parser.add_argument(
        "-p",
        "--topology",
        required=True,
        help="Topology/PDB file",
    )

    parser.add_argument(
        "-t",
        "--trajectory",
        required=True,
        help="Trajectory file",
    )

    parser.add_argument(
        "--base",
        default="A",
        help="Base bead atom name (default: A)",
    )

    parser.add_argument(
        "--start",
        type=int,
        default=1000,
        help="First trajectory frame (default: 1000)",
    )

    parser.add_argument(
        "--stop",
        type=int,
        default=None,
        help="Final trajectory frame",
    )

    parser.add_argument(
        "--step",
        type=int,
        default=1,
        help="Trajectory stride (default: 1)",
    )

    parser.add_argument(
        "-o",
        "--output",
        default="stack_propensity.dat",
        help="Output file",
    )

    # optional wrapping

    parser.add_argument(
        "--wrap",
        action="store_true",
        help=(
            "Apply unwrap -> center RNA -> wrap "
            "trajectory transformations"
        ),
    )

    parser.add_argument(
        "--box-size",
        type=float,
        default=None,
        help=(
            "Cubic box length in Angstrom. "
            "Required when --wrap is used."
        ),
    )

    parser.add_argument(
        "--rna-selection",
        type=str,
        default=None,
        help=(
            "MDAnalysis selection for RNA centering, "
            'e.g. "resname ADE"'
        ),
    )

    args = parser.parse_args()

    analyze_stacking(
        topology=args.topology,
        trajectory=args.trajectory,
        start_frame=args.start,
        stop_frame=args.stop,
        step=args.step,
        nucleotide_name=args.base,
        output=args.output,
        apply_wrap=args.wrap,
        box_size=args.box_size,
        rna_selection=args.rna_selection,
    )


if __name__ == "__main__":
    main()


# Examples
# Without wrapping
# python stacking_analysis.py \
#     -p rA30-Mg1-Na20-1.pdb \
#     -t md3.dcd \
#     --base A \
#     --start 1000 \
#     --step 10


# With wrapping
# python stacking_analysis.py \
#     -p rA30-Mg1-Na20-1.pdb \
#     -t md3.dcd \
#     --base A \
#     --start 1000 \
#     --step 10 \
#     --wrap \
#     --box-size 720 \
#     --rna-selection "resname ADE"
