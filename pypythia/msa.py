import math
import os
import random
import statistics
from collections import Counter
from itertools import product
from tempfile import NamedTemporaryFile
import warnings

from Bio import AlignIO, SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import DistanceCalculator

from pypythia.custom_types import *
from pypythia.custom_errors import PyPythiaException

STATE_CHARS = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ!\"#$%&'()*+,/:;<=>@[\\]^_{|}~"
DNA_CHARS = "ATUCGMRWSYKVHDBN"
GAP_CHARS = "-?."


class MSA:
    """Class structure for Multiple Sequence Alignemnt statistics.

    This class provides methods for computing MSA attributes.

    Args:
        msa_file (FilePath): Path to the MSA file the statistics should be computed for. The file must be either in "fasta" or "phylip" format.
        msa_name (str, optional): Better readable name for the MSA. If not set, the name of the file is used.
        file_format (FileFormat, optional): If not set, Pythia attempts to autodetect the file format.

    Attributes:
        msa_file (FilePath): Path to the corresponding MSA file.
        msa_name (str): Name of the MSA. Can be either set or is inferred automatically based on the msa_file.
        data_type (DataType): Data type of the MSA.
        file_format (FileFormat): File format of the MSA.
        msa (MultipleSeqAlignment): Biopython MultipleSeqAlignment object for the given msa_file.

    Raises:
        ValueError: If the file format of the given MSA is not FASTA or PHYLIP.
        ValueError: If the data type of the given MSA cannot be inferred.
    """

    def __init__(self, msa_file: FilePath, msa_name: str = None, file_format: FileFormat = None):
        self.msa_file = msa_file
        self.file_format = file_format if file_format is not None else self._get_file_format()
        self.data_type = self.guess_data_type()

        if msa_name:
            self.msa_name = msa_name
        else:
            self.msa_name = os.path.split(msa_file)[1]

        self._set_msa_object(self.msa_file)

    def _set_msa_object(self, msa_file: FilePath):
        with NamedTemporaryFile(mode="w") as tmpfile:
            self._convert_msa_to_biopython_format(msa_file, tmpfile)
            tmpfile.flush()
            try:
                self.msa = AlignIO.read(tmpfile.name, format=self.file_format.value)
            except Exception as e:
                raise PyPythiaException("Error reading the provided MSA: ", msa_file) from e

    def _get_updated_sequence(self, sequence: str, char_mapping: dict):
        sequence = sequence.upper()
        for char, new_char in char_mapping.items():
            sequence = sequence.replace(char, new_char)

        return sequence

    def _replace_phylip_sequence_chars(self, msa_file: FilePath, char_mapping: dict, new_file: FilePath):
        msa_content = open(msa_file)
        new_content = []

        # phylip file contains the number of taxa and sites as integer numbers in the first line
        first_line = msa_content.readline().strip()
        _ntaxa, _nsites, *_ = first_line.split()
        new_content.append(first_line + "\n")

        for i, line in enumerate(msa_content):
            if i < int(_ntaxa):
                # for the first _ntaxa lines: replace all chars only in the sequences
                taxon, sequence = line.split(maxsplit=1)
                sequence = self._get_updated_sequence(sequence, char_mapping)
                new_line = f"{taxon} {sequence}"
            else:
                new_line = self._get_updated_sequence(line, char_mapping)

            new_content.append(new_line)

        msa_content.close()

        new_file.write("".join(new_content))

    def _replace_fasta_sequence_chars(self, msa_file: FilePath, char_mapping: dict, new_file: FilePath):
        msa_content = open(msa_file)
        new_content = []

        for line in msa_content:
            if line.startswith(">"):
                # this line contains the taxon name, we can simply add it
                new_line = line
            else:
                new_line = self._get_updated_sequence(line, char_mapping)

            new_content.append(new_line)

        msa_content.close()

        new_file.write("".join(new_content))

    def _convert_msa_to_biopython_format(self, msa_file: FilePath, tmpfile: NamedTemporaryFile) -> None:
        """
        The unknown char in DNA MSA files for Biopython to work
        has to be "-" instead of "X" or "N" -> replace all occurrences
        All "?" are gaps -> convert to "-"
        Also all "U" need to be "T"

        For AA and MORPH files, we only have to replace the "?" gap symbol with "-"

        We do however, need to be careful as to not replace any of these chars in the taxon names.
        """
        if self.data_type == DataType.DNA:
            char_mapping = {
                "X": "-",
                "N": "-",
                "?": "-",
                "U": "T"
            }
        else:
            char_mapping = {
                "?": "-"
            }

        if self.file_format == FileFormat.PHYLIP:
            self._replace_phylip_sequence_chars(msa_file, char_mapping, tmpfile)
        elif self.file_format == FileFormat.FASTA:
            self._replace_fasta_sequence_chars(msa_file, char_mapping, tmpfile)
        else:
            raise PyPythiaException("Unsupported file format: ", self.file_format)

    def _get_file_format(self) -> FileFormat:
        first_line = open(self.msa_file).readline().strip()

        if first_line.startswith(">"):
            return FileFormat.FASTA

        try:
            # phylip file contains two integer numbers in the first line separated by whitespace
            _num1, _num2, *_ = first_line.split()
            _num1 = int(_num1)
            _num2 = int(_num2)
            # in case these conversions worked, the file is (most likely) in phylip format
            return FileFormat.PHYLIP
        except:
            raise ValueError(
                f"The file type of {self.msa_file} could not be autodetected, please check that the file contains data in phylip or fasta format."
            )

    def _get_unique_sequences(self):
        unique_sequences = set()
        unique_seq_objects = []

        for seq_object in self.msa:
            sequence = seq_object.seq

            if sequence not in unique_sequences:
                unique_sequences.add(sequence)
                unique_seq_objects.append(seq_object)

        return unique_seq_objects

    def save_reduced_alignment(self, reduced_msa_file: FilePath, replace_original: bool = False):
        """ Removes all duplicate sequences from the MSA and stores the reduced alignment in reduced_msa_file.
        In case of duplicate sequences, the reduced alignment only contains the first occurrence of the respective sequence.

        Args:
            reduced_msa_file (FilePath): Filename where to store the reduced alignment.
            replace_original (bool): Optional switch.
                If True, self.msa_file is replaced with the reduced alignment and
                self.msa with the MultipleSeqAlignment object of the reduced MSA.
                Self.msa_name will be appended with "_reduced"
        """
        contains_duplicates = self.contains_duplicate_sequences()

        if not contains_duplicates:
            raise PyPythiaException("This MSA does not contain duplicate sequences, MSA reduction does not make sense.")

        unique_sequence_objects = self._get_unique_sequences()
        SeqIO.write(unique_sequence_objects, reduced_msa_file, self.file_format.value)

        if replace_original:
            self._set_msa_object(reduced_msa_file)
            self.msa_name += "_reduced"
            self.msa_file = reduced_msa_file

    def contains_duplicate_sequences(self) -> bool:
        num_unique = len({s.seq for s in self.msa})

        return self.number_of_taxa() > num_unique

    def guess_data_type(self) -> DataType:
        """Guesses and returns the data type (DNA, AA, Morphological data) based on the contents of the MSA.
        The data type is guessed as DNA, if the sequences only contains DNA_CHARS or GAP_CHARS.
        The data type is guesses as morphological if the sequences contain digits.

        Returns:
            data_type (DataType): A guess of the data type for this MSA.
        """
        format = self.file_format
        msa_content = open(self.msa_file).readlines()
        sequence_chars = set()

        if format == FileFormat.PHYLIP:
            # the sequences start in the second line and the schema is "taxon_name [SEQUENCE]"
            for line in msa_content[1:]:
                line = line.strip()
                if not line:
                    continue
                try:
                    # ...with some files the sequences are split into multiple lines
                    # in this case the taxon name is skipped and splitting will raise a ValueError
                    # if this happens, we simply add the rest of the block to the set
                    _, sequence = line.split(None, 1)
                except ValueError:
                    # read the line as is
                    sequence = line
                # remove whitespace and add the characters to the set of unique characters of this MSA
                sequence = sequence.strip().replace(" ", "")
                for char in set(sequence):
                    sequence_chars.add(char.upper())
        elif format == FileFormat.FASTA:
            # taxon names start with a ">", treat every other line as sequence line
            for line in msa_content:
                line = line.strip()
                if line.startswith(">"):
                    continue
                else:
                    for char in set(line):
                        sequence_chars.add(char.upper())
        else:
            raise ValueError(
                f"Unsupported MSA file format {format}. Supported formats are phylip and fasta."
            )

        is_dna = all([(c in DNA_CHARS) or (c in GAP_CHARS) for c in sequence_chars])
        is_morph = any([c.isdigit() for c in sequence_chars])

        # now check whether the sequence_chars contain only DNA and GAP chars or not
        if is_dna:
            return DataType.DNA
        elif is_morph:
            return DataType.MORPH
        else:
            return DataType.AA

    def get_raxmlng_model(self) -> Model:
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
            unique_states = set()
            for seq in self.msa:
                seq = str(seq.seq).replace("-", "")
                unique_states = unique_states.union(set(seq))
            # the number of unique states is irrelevant for RAxML-NG, it only cares about the max state value...
            num_states = int(max(unique_states)) + 1
            return f"MULTI{num_states}_GTR"
        else:
            raise ValueError("Unsupported data type: ", self.data_type)

    def number_of_taxa(self) -> int:
        """Returns the number of taxa of the MSA.

        Returns:
            n_taxa (int): The number of taxa of the MSA
        """
        return len(self.msa)

    def number_of_sites(self) -> int:
        """Returns the number of sites of the MSA.

        Returns:
            n_sites (int): The number of sites of the MSA
        """
        return self.msa.get_alignment_length()

    def column_entropies(self) -> List[float]:
        """Returns the shannon entropy (in bits) for each site in the MSA.

        Returns:
            column_entropies (List[float]): List of shannon entropies corresponding to the MSA site. Each per-site entropy is >= 0.
        """

        def _remove_gaps_from_sequence(seq):
            for char in GAP_CHARS:
                seq = seq.replace(char, "")
            return seq

        entropies = []
        for i in range(self.msa.get_alignment_length()):
            column = _remove_gaps_from_sequence(self.msa[:, i]).upper()
            entropy = 0

            for char in STATE_CHARS:
                count = str.count(column, char)
                if count == 0:
                    entropy_x = 0
                else:
                    prob = count / len(column)
                    entropy_x = prob * math.log2(prob)

                entropy += entropy_x

            entropy = -entropy

            assert (
                    entropy >= 0
            ), f"Entropy negative, check computation. Entropy is {entropy}"

            entropies.append(entropy)
        return entropies

    def entropy(self) -> float:
        """Returns the shannon entropy (in bits) of the MSA.

        Returns:
            entropy (float): Shannon entropy of the MSA. The entropy is >= 0.
        """
        return statistics.mean(self.column_entropies())

    def pattern_entropy(self) -> float:
        """Returns an entropy like metric based on the number and frequency of patterns in the MSA.

        Returns:
            pattern entropy (float): The pattern entropy like measure. The pattern entropy is >= 0.
        TODO: Documentation and Test
        """
        msa_length = self.number_of_sites()

        sites = []
        for i in range(msa_length):
            sites.append(self.msa[:, i])

        site_counts = Counter(sites)
        mult = 0
        for N_i in site_counts.values():
            mult += N_i * math.log(N_i)

        return mult

    def bollback_multinomial(self) -> float:
        """Returns the bollback multinomial statistic of the MSA.

        According to Bollback, JP: Bayesian model adequacy and choice in phylogenetics (2002)

        Returns:
            bollback (float): The bollback multionomial statistic of the MSA. The bollback multinomial statistic is <= 0.
        """
        pattern_entropy = self.pattern_entropy()
        msa_length = self.number_of_sites()
        bollback = pattern_entropy - msa_length * math.log(msa_length)
        return bollback

    def _get_distance_matrix(self, num_samples: int):
        """
        For large MSAs (i.e. more than num_samples taxa), computing the distance matrix
        is computationally very expensive.
        So for large MSAs, we rather compute the distance matrix on a subsample of at most num_samples sequences
        """
        if self.number_of_taxa() > num_samples:
            sample_population = range(self.number_of_taxa())
            selection = sorted(random.sample(sample_population, num_samples))
            _msa = MultipleSeqAlignment([self.msa[el] for el in selection])
        else:
            _msa = self.msa

        model = "blastn" if self.data_type == DataType.DNA else "blosum62"
        calculator = DistanceCalculator(model=model)
        return calculator.get_distance(_msa)

    def treelikeness_score(self, num_samples: int = 100) -> float:
        """Returns the treelikeness score of the MSA.

        Computed according to
        δ Plots: A Tool for Analyzing Phylogenetic Distance Data, Holland, Huber, Dress and Moulton (2002)
        https://doi.org/10.1093/oxfordjournals.molbev.a004030

        For large MSAs (i.e. more than num_samples taxa), computing the distance matrix
        is computationally very expensive.
        So for large MSAs, we rather compute the distance matrix on a subsample of at most num_samples sequences.

        Args:
            num_samples (int): When computing the distance matrix required for this statistic compute at most min(n_taxa, num_samples) pairwise distances.

        Returns:
            treelikeness (float): Treelikeness of the MSA. The treelikeness is in the value range [0.0, 1.0].
                The lower the treelikeness the stronger the phylogenetic signal of the MSA.
        """
        if self.data_type == DataType.MORPH:
            warnings.warn("Computing the treelikeness score for morphological data is currently not supported.")
            return -1.0

        num_samples = min(self.number_of_taxa(), num_samples)
        dm = self._get_distance_matrix(num_samples)

        options = list(range(len(dm)))

        frac = num_samples // 4
        X = options[:frac]
        Y = options[frac: 2 * frac]
        U = options[2 * frac: 3 * frac]
        V = options[3 * frac:]

        res = product(X, Y, U, V)
        deltas = []
        for x, y, u, v in res:
            dxv = abs(dm[x, v])
            dyu = abs(dm[y, u])
            dxu = abs(dm[x, u])
            dyv = abs(dm[y, v])
            dxy = abs(dm[x, y])
            duv = abs(dm[u, v])

            dxv_yu = dxv + dyu
            dxu_yv = dxu + dyv
            dxy_uv = dxy + duv

            vals = sorted([dxv_yu, dxu_yv, dxy_uv])
            smallest = vals[0]
            intermediate = vals[1]
            largest = vals[2]

            numerator = largest - intermediate
            denominator = largest - smallest

            if denominator == 0:
                delta = 0
            else:
                delta = numerator / denominator
            assert delta >= 0
            assert delta <= 1
            deltas.append(delta)

        return statistics.mean(deltas)
