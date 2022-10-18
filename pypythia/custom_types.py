from typing import List, Tuple, Dict, Union
from enum import Enum

FilePath = str
Command = List[str]
Model = Union[str, FilePath]
Executable = str


class DataType(Enum):
    """Data type for MSAs.
    - DNA = DNA data
    - AA = Protein data
    - MORPH = morphological data
    """
    DNA = "DNA"
    AA = "AA"
    MORPH = "MORPH"


class FileFormat(Enum):
    """File formats for MSAs
    - PHYLIP = phylip-relaxed
    - FASTA = fasta
    """
    PHYLIP = "phylip-relaxed"
    FASTA = "fasta"