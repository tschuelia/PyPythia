from custom_types import *
import subprocess


def run_cmd(cmd: Command) -> None:
    try:
        subprocess.check_output(cmd)
    except Exception as e:
        print(f"Error running command \"{' '.join(cmd)}\"")
        raise e


def get_value_from_line(line: str, search_string: str) -> float:
    line = line.strip()
    if search_string in line:
        _, value = line.rsplit(" ", 1)
        return float(value)

    raise ValueError(
        f"The given line '{line}' does not contain the search string '{search_string}'."
    )


def get_single_value_from_file(input_file: FilePath, search_string: str) -> float:
    with open(input_file) as f:
        lines = f.readlines()

    for l in lines:
        if search_string in l:
            return get_value_from_line(l, search_string)

    raise ValueError(
        f"The given input file {input_file} does not contain the search string '{search_string}'."
    )


def get_raxmlng_abs_rf_distance(log_file: FilePath) -> float:
    STR = "Average absolute RF distance in this tree set:"
    return get_single_value_from_file(log_file, STR)


def get_raxmlng_rel_rf_distance(log_file: FilePath) -> float:
    STR = "Average relative RF distance in this tree set:"
    return get_single_value_from_file(log_file, STR)


def get_raxmlng_num_unique_topos(log_file: FilePath) -> int:
    STR = "Number of unique topologies in this tree set:"
    return int(get_single_value_from_file(log_file, STR))


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
