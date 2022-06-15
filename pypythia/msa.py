import math
import os
import random
import statistics
from collections import Counter
from itertools import product
from tempfile import NamedTemporaryFile

from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import DistanceCalculator

from pypythia.custom_types import *

STATE_CHARS = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ!\"#$%&'()*+,/:;<=>@[\\]^_{|}~"
DNA_CHARS = "ATUCGMRWSYKVHDBN"
GAP_CHARS = "-?."


class MSA:
    """Class structure for Multiple Sequence Alignemnt statistics.

    This class provides methods for computing MSA attributes.

    Args:
        msa_file (str): Path to the MSA file the statistics should be computed for. The file must be either in "fasta" or "phylip" format.

    Attributes:
        msa_file (str): Path to the corresponding MSA file.
        msa_name (str): Name of the MSA. Can be either set or is inferred automatically based on the msa_file.
        data_type (str): Data type of the MSA. Can be either "DNA" for DNA data or "AA" for protein data.
        msa (MultipleSeqAlignment): Biopython MultipleSeqAlignment object for the given msa_file.

    Raises:
        ValueError: If the file format of the given MSA is not "fasta" or "phylip".
        ValueError: If the data type of the given MSA cannot be inferred.
    """

    def __init__(self, msa_file: FilePath, msa_name: str = None):
        self.msa_file = msa_file
        self.data_type = self.guess_data_type()

        if msa_name:
            self.msa_name = msa_name
        else:
            self.msa_name = os.path.split(msa_file)[1]

        with NamedTemporaryFile(mode="w") as tmpfile:
            if self.data_type == "DNA":
                self._convert_dna_msa_to_biopython_format(tmpfile)
            else:
                self._convert_aa_msa_to_biopython_format(tmpfile)

            tmpfile.flush()
            self.msa = AlignIO.read(tmpfile.name, format=self._get_file_format())

    def _convert_dna_msa_to_biopython_format(self, tmpfile: NamedTemporaryFile) -> None:
        """
        The unknonwn char in DNA MSA files for Biopython to work
        has to be "-" instead of "X" or "N" -> replace all occurences
        All "?" are gaps -> convert to "-"
        Also all "U" need to be "T"
        """
        with open(self.msa_file) as f:
            repl = f.read().replace("X", "-")
            repl = repl.replace("N", "-")
            repl = repl.replace("?", "-")
            repl = repl.replace("U", "T")
            repl = repl.upper()

        tmpfile.write(repl)

    def _convert_aa_msa_to_biopython_format(self, tmpfile: NamedTemporaryFile) -> None:
        """
        The unknonwn char in AA MSA files for Biopython to work
        All "?" are gaps -> convert to "-"
        """
        with open(self.msa_file) as f:
            repl = f.read().replace("?", "-")
            repl = repl.upper()

        tmpfile.write(repl)

    def _get_file_format(self) -> FileFormat:
        first_line = open(self.msa_file).readline().strip()

        try:
            # phylip file contains two integer numbers in the first line separated by whitespace
            _num1, _num2, *_ = first_line.split()
            _num1 = int(_num1)
            _num2 = int(_num2)
            # in case these conversions worked, the file is (most likely) in phylip format
            return "phylip-relaxed"
        except:
            # if the MSA is in fasta format, the first line should start with a '>' character
            if first_line.startswith(">"):
                return "fasta"

        raise ValueError(
            f"The file type of this MSA could not be autodetected, please check file."
        )

    def guess_data_type(self) -> DataType:
        format = self._get_file_format()
        msa_content = open(self.msa_file).readlines()

        sequence_chars = set()

        if format == "phylip-relaxed":
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
                    sequence_chars.add(char)
        elif format == "fasta":
            seen_taxon_name = False
            # taxon names start with a ">" followed by the sequence
            for line in msa_content:
                line = line.strip()
                if line.startswith(">"):
                    seen_taxon_name = True
                    continue
                else:
                    if seen_taxon_name:
                        for char in set(line):
                            sequence_chars.add(char)
                        seen_taxon_name = False
        else:
            raise ValueError(
                f"Unsupported MSA file format {format}. Supported formats are phylip and fasta."
            )

        # now check whether the sequence_chars contain only DNA and GAP chars or not
        if all([(c in DNA_CHARS) or (c in GAP_CHARS) for c in sequence_chars]):
            return "DNA"
        else:
            return "AA"

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

    def bollback_multinomial(self) -> float:
        """Returns the bollback multinomial statistic of the MSA.

        According to Bollback, JP: Bayesian model adequacy and choice in phylogenetics (2002)

        Returns:
            bollback (float): The bollback multionomial statistic of the MSA. The bollback multinomial statistic is <= 0.
        """
        msa_length = self.number_of_sites()

        sites = []
        for i in range(msa_length):
            sites.append(self.msa[:, i])

        site_counts = Counter(sites)
        mult = 0
        for i in site_counts:
            N_i = site_counts[i]
            mult += N_i * math.log(N_i)

        mult = mult - msa_length * math.log(msa_length)
        return mult

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

        model = "blastn" if self.data_type == "DNA" else "blosum62"
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
