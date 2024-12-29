import sys
import textwrap
import time

import loguru

from pypythia import __version__

SCRIPT_START = time.perf_counter()


logger = loguru.logger
logger.remove()
logger.add(sys.stderr, format="{message}")


def get_header():
    return textwrap.dedent(
        f"PyPythia version {__version__} released by The Exelixis Lab\n"
        f"Developed by: Julia Haag\n"
        f"Latest version: https://github.com/tschuelia/PyPythia\n"
        f"Questions/problems/suggestions? Please open an issue on GitHub.\n",
    )


def log_runtime_information(message, log_runtime=True):
    if log_runtime:
        seconds = time.perf_counter() - SCRIPT_START
        fmt_time = time.strftime("%H:%M:%S", time.gmtime(seconds))
        time_string = f"[{fmt_time}] "
    else:
        time_string = ""
    logger.info(f"{time_string}{message}")
