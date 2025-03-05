import importlib.metadata

try:
    __version__ = importlib.metadata.distribution(__name__).version
except Exception:
    __version__ = "unknown"

# Required if the package was installed via PyPi...
if __version__ == "unknown":
    try:
        from importlib.metadata import version

        __version__ = version("PythiaPhyloPredictor")
    except Exception:
        pass
