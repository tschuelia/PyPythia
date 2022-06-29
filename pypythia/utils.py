from pypythia.custom_types import *
import subprocess


def run_cmd(cmd: Command) -> None:
    try:
        subprocess.check_output(cmd)
    except Exception as e:
        raise RuntimeError(f"Error running command \"{' '.join(cmd)}\"") from e


def get_value_from_line(line: str, search_string: str) -> float:
    line = line.strip()
    if search_string in line:
        _, value = line.rsplit(" ", 1)
        return float(value)

    raise ValueError(
        f"The given line '{line}' does not contain the search string '{search_string}'."
    )

