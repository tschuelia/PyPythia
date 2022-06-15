from pypythia.custom_types import *
from pypythia.utils import get_value_from_line


def get_patterns_gaps_invariant(log_file: FilePath) -> Tuple[int, float, float]:
    """Method that parses the number of patterns, proportion of gaps, and proportion of invariant sites in the given log_file.

    Args:
        log_file (str): Filepath of a RAxML-NG log file.

    Returns:
        n_patterns (int): Number of unique patterns in the given MSA.
        prop_gaps (float): Proportion of gaps in the given MSA.
        prop_inv (float): Proportion of invariant sites in the given MSA.

    Raises:
        ValueError: If the given log file does not contain the number of patterns, proportion of gaps or proportion of invariant sites.
    """
    patterns = None
    gaps = None
    invariant = None
    for line in open(log_file).readlines():
        if line.startswith("Alignment sites"):
            # number of alignment patterns
            # Alignment sites / patterns: 1940 / 933
            _, numbers = line.split(":")
            _, patterns = [int(el) for el in numbers.split("/")]
        elif line.startswith("Gaps"):
            # proportion of gaps
            _, number = line.split(":")
            percentage, _ = number.strip().split(" ")
            gaps = float(percentage) / 100.0
        elif line.startswith("Invariant sites"):
            # proportion invariant sites
            _, number = line.split(":")
            percentage, _ = number.strip().split(" ")
            invariant = float(percentage) / 100.0

    if patterns is None or gaps is None or invariant is None:
        raise ValueError("Error parsing raxml-ng log")

    return patterns, gaps, invariant


def get_raxmlng_rfdist_results(log_file: FilePath) -> Tuple[float, float, float]:
    """Method that parses the number of unique topologies, relative RF-Distance, and absolute RF-Distance in the given log file.

    Args:
        log_file (str): Filepath of a RAxML-NG log file.

    Returns:
        num_topos (int): Number of unique topologies of the given set of trees.
        rel_rfdist (float): Relative RF-Distance of the given set of trees. Computed as average over all pairwise RF-Distances. Value between 0.0 and 1.0.
        abs_rfdist (float): Absolute RF-Distance of the given set of trees.

    Raises:
        ValueError: If the given log file does not contain the unique topologies, relative RF-Distance, or absolute RF-Distance.
    """
    abs_rfdist = None
    rel_rfdist = None
    num_topos = None

    for line in open(log_file).readlines():
        line = line.strip()

        if "Average absolute RF distance in this tree set:" in line:
            abs_rfdist = get_value_from_line(
                line, "Average absolute RF distance in this tree set:"
            )
        elif "Average relative RF distance in this tree set:" in line:
            rel_rfdist = get_value_from_line(
                line, "Average relative RF distance in this tree set:"
            )
        elif "Number of unique topologies in this tree set:" in line:
            num_topos = get_value_from_line(
                line, "Number of unique topologies in this tree set:"
            )

    if abs_rfdist is None or rel_rfdist is None or num_topos is None:
        raise ValueError("Error parsing raxml-ng log.")

    return num_topos, rel_rfdist, abs_rfdist
