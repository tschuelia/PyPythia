import math
import pathlib
from collections import Counter
from io import StringIO
from typing import Optional

import numpy as np
import numpy.typing as npt
from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from pypythia.custom_errors import PyPythiaException
from pypythia.custom_types import DataType, FileFormat

GAP = b"-"
GAP_ORD = ord(GAP)
GAP_CHARS = [GAP, b"?", b".", b"X", b"*"]

NUCLEOTIDES = [b"A", b"C", b"G", b"T"]
DNA_GAP_CHARS = GAP_CHARS + [b"N"]

DNA_AMBIGUITY_CODES = {
    b"R": [b"A", b"G"],
    b"K": [b"G", b"T"],
    b"S": [b"G", b"C"],
    b"Y": [b"C", b"T"],
    b"M": [b"A", b"C"],
    b"W": [b"A", b"T"],
    b"B": [b"C", b"G", b"T"],
    b"H": [b"A", b"C", b"T"],
    b"D": [b"A", b"G", b"T"],
    b"V": [b"A", b"C", b"G"],
}

DNA_AMBIGUITY_CHARS = list(DNA_AMBIGUITY_CODES.keys())

DNA_AMBIGUITY_MAP = {
    nt: set(
        map(
            ord,
            [GAP, nt] + [c for c, vals in DNA_AMBIGUITY_CODES.items() if nt in vals],
        )
    )
    for nt in NUCLEOTIDES
}

DNA_CHARS = NUCLEOTIDES + list(DNA_AMBIGUITY_CODES.keys()) + DNA_GAP_CHARS


AMINO_ACIDS = [
    b"A",
    b"C",
    b"D",
    b"E",
    b"F",
    b"G",
    b"H",
    b"I",
    b"K",
    b"L",
    b"M",
    b"N",
    b"P",
    b"Q",
    b"R",
    b"S",
    b"T",
    b"V",
    b"W",
    b"Y",
]
AA_AMBIGUITY_CODES = {
    b"B": [b"D", b"N"],
    b"Z": [b"E", b"Q"],
}

AA_AMBIGUITY_CHARS = list(AA_AMBIGUITY_CODES.keys())

AA_AMBIGUITY_MAP = {
    aa: set(
        map(
            ord, [GAP, aa] + [c for c, vals in AA_AMBIGUITY_CODES.items() if aa in vals]
        )
    )
    for aa in AMINO_ACIDS
}

AA_CHARS = AMINO_ACIDS + list(AA_AMBIGUITY_CODES.keys()) + GAP_CHARS


def _get_file_format(msa_file: pathlib.Path) -> FileFormat:
    first_line = msa_file.open().readline().strip()

    if first_line.startswith(">"):
        return FileFormat.FASTA

    try:
        # phylip file contains two integer numbers in the first line separated by whitespace
        _num1, _num2, *_ = first_line.split()
        _num1 = int(_num1)
        _num2 = int(_num2)
        # in case these conversions worked, the file is (most likely) in phylip format
        return FileFormat.PHYLIP
    except Exception:
        raise ValueError(
            f"The file type of {msa_file} could not be determined. "
            f"Please check that the file contains data in phylip or fasta format."
        )


def _guess_dtype(sequences: npt.NDArray):
    seq_chars = np.unique(sequences)

    # First, check if any character is a digit, if yes: morphological data
    if np.any(np.char.isdigit(seq_chars)):
        return DataType.MORPH

    # Next, check if all characters are valid DNA_CHARS, if yes: DNA data
    if np.all(np.isin(seq_chars, DNA_CHARS)):
        return DataType.DNA

    # Next, check if all characters are valid AA_CHARS, if yes: AA data
    if np.all(np.isin(seq_chars, AA_CHARS)):
        return DataType.AA

    raise ValueError(
        f"Data type for character set could not be inferred: {[c.decode() for c in seq_chars]}."
        f" Invalid characters for DNA: {[c.decode() for c in set(seq_chars) - set(DNA_CHARS)]}."
        f" Invalid characters for AA: {[c.decode() for c in set(seq_chars) - set(AA_CHARS)]}."
    )


class MSA:
    def __init__(
        self, taxa: npt.NDArray, sequences: npt.NDArray, data_type: DataType, name: str
    ):
        if taxa.shape[0] != sequences.shape[0]:
            raise PyPythiaException(
                "Number of taxa and sequences do not match: "
                f"{taxa.shape[0]} != {sequences.shape[0]}."
            )

        self.taxa = taxa
        self.sequences = sequences
        self.data_type = data_type
        self.name = name

        self.n_taxa, self.n_sites = self.sequences.shape

    def __str__(self):
        return f"MSA(name={self.name}, n_taxa={self.n_taxa}, n_sites={self.n_sites}, data_type={self.data_type})"

    def __repr__(self):
        return str(self)

    def contains_full_gap_sequences(self):
        return np.any(np.all(self.sequences == GAP, axis=1))

    def contains_duplicate_sequences(self):
        unique_sequences = np.unique(self.sequences, axis=0)
        return unique_sequences.shape[0] < self.sequences.shape[0]

    @property
    def n_patterns(self):
        un = set([c.tobytes() for c in self.sequences.T])
        return len(un) - ((GAP * self.n_taxa) in un)

    @property
    def percentage_gaps(self):
        full_gap_sites_removed = self.sequences[
            :, ~np.all(self.sequences.T == GAP, axis=1)
        ]
        return np.sum(full_gap_sites_removed == GAP) / full_gap_sites_removed.size

    @property
    def percentage_invariant(self):
        if self.data_type == DataType.DNA:
            charmap = DNA_AMBIGUITY_MAP
        elif self.data_type == DataType.AA:
            charmap = AA_AMBIGUITY_MAP
        else:
            charmap = {}

        non_gap_site_count = 0
        invariant_count = 0

        for site in self.sequences.T:
            site = set(site.tobytes())
            if site == {GAP_ORD}:
                # full-gap sites are not counted as invariant
                continue

            non_gap_site_count += 1

            for allowed in charmap.values():
                if site.issubset(allowed):
                    invariant_count += 1
                    break

        return invariant_count / non_gap_site_count

    def entropy(self):
        def _site_entropy(site):
            site_counter = Counter(site.tobytes())
            site_counter.pop(GAP_ORD, None)

            counts = np.array(list(site_counter.values()))
            probabilities = counts / np.sum(counts)
            return -np.sum(probabilities * np.log2(probabilities))

        return np.mean([_site_entropy(site) for site in self.sequences.T])

    def pattern_entropy(self):
        patterns = [c.tobytes() for c in self.sequences.T]
        pattern_counter = Counter(patterns)
        pattern_counts = np.array(list(pattern_counter.values()))
        return np.sum(pattern_counts * np.log(pattern_counts))

    def bollback_multinomial(self):
        pattern_entropy = self.pattern_entropy()
        return pattern_entropy - self.n_sites * math.log(self.n_sites)

    def get_raxmlng_model(self) -> str:
        """Returns a RAxML-NG model string based on the data type

        Returns:
             model_string (string): RAxML-NG model string
                For DNA data: GTR+G
                For Protein (AA) data: LG+G
                For morphological data: MULTIx_GTR where x refers to the maximum state value in the alignment
        """
        if self.data_type == DataType.DNA:
            return "GTR+G"
        elif self.data_type == DataType.AA:
            return "LG+G"
        elif self.data_type == DataType.MORPH:
            unique = np.unique(self.sequences)
            # the number of unique states is irrelevant for RAxML-NG, it only cares about the max state value...
            num_states = int(max(unique)) + 1
            return f"MULTI{num_states}_GTR"
        else:
            raise ValueError("Unsupported data type: ", self.data_type)

    def save(
        self, output_file: pathlib.Path, file_format: FileFormat = FileFormat.PHYLIP
    ):
        _biopython_sequences = [
            SeqRecord(Seq(seq.tobytes()), id=taxon)
            for seq, taxon in zip(self.sequences, self.taxa)
        ]
        SeqIO.write(_biopython_sequences, output_file, file_format.value)


def parse(
    msa_file: pathlib.Path,
    file_format: Optional[FileFormat] = None,
    data_type: Optional[DataType] = None,
) -> MSA:
    file_format = file_format or _get_file_format(msa_file)
    msa_content = StringIO(msa_file.read_text().upper())
    _msa = AlignIO.read(msa_content, format=file_format.value)
    sequences = (
        np.frombuffer(b"".join([rec.seq._data for rec in _msa]), dtype="S1")
        .reshape(len(_msa), -1)
        .copy()
    )
    taxon_names = np.array([rec.id for rec in _msa])

    if not data_type:
        data_type = _guess_dtype(sequences)

    char_mapping = {c: GAP for c in GAP_CHARS}
    if data_type == DataType.DNA:
        # replace all U characters with a T for convenience
        char_mapping.update({b"U": b"T"})
        char_mapping.update({c: GAP for c in DNA_GAP_CHARS})

    for old_char, new_char in char_mapping.items():
        sequences[sequences == old_char] = new_char

    return MSA(taxon_names, sequences, data_type, msa_file.name)


def remove_full_gap_sequences(msa: MSA, msa_name: Optional[str] = None) -> MSA:
    if not msa.contains_full_gap_sequences():
        raise PyPythiaException("No full-gap sequences found in MSA.")

    is_full_gap_sequence = np.all(msa.sequences == GAP, axis=1)
    non_full_gap_sequences = msa.sequences[~is_full_gap_sequence]
    non_full_gap_taxa = msa.taxa[~is_full_gap_sequence]

    return MSA(
        non_full_gap_taxa, non_full_gap_sequences, msa.data_type, msa_name or msa.name
    )


def deduplicate_sequences(msa: MSA, msa_name: Optional[str] = None) -> MSA:
    if not msa.contains_duplicate_sequences():
        raise PyPythiaException("No duplicates found in MSA.")

    unique_sequences, unique_indices = np.unique(
        msa.sequences, axis=0, return_index=True
    )
    unique_taxa = msa.taxa[unique_indices]

    return MSA(unique_taxa, unique_sequences, msa.data_type, msa_name or msa.name)
