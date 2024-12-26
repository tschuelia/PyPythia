import importlib.metadata

try:
    __version__ = importlib.metadata.distribution(__name__).version
except Exception:
    __version__ = "unknown"
