from custom_types import *
from utils import get_single_value_from_file


def get_patterns_gaps_invariant(log_file: FilePath) -> Tuple[int, float, float]:
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


def get_raxmlng_abs_rf_distance(log_file: FilePath) -> float:
    STR = "Average absolute RF distance in this tree set:"
    return get_single_value_from_file(log_file, STR)


def get_raxmlng_rel_rf_distance(log_file: FilePath) -> float:
    STR = "Average relative RF distance in this tree set:"
    return get_single_value_from_file(log_file, STR)


def get_raxmlng_num_unique_topos(log_file: FilePath) -> int:
    STR = "Number of unique topologies in this tree set:"
    return int(get_single_value_from_file(log_file, STR))
