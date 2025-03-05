import pathlib
import shutil

DEFAULT_RAXMLNG_EXE = (
    pathlib.Path(shutil.which("raxml-ng")) if shutil.which("raxml-ng") else None
)
DEFAULT_MODEL_FILE = pathlib.Path(__file__).parent / "predictors" / "latest.txt"
