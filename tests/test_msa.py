from tests.fixtures import *
from pypythia.custom_types import DataType, FileFormat
from pypythia.msa import MSA
from pypythia.custom_errors import PyPythiaException


class TestMSAFeatures:
    def test_contains_duplicate_sequences(self, msa_with_duplicate_sequences, msa_without_duplicate_sequences):
        # small_msa_with_signal does not contain duplicates
        assert msa_with_duplicate_sequences.contains_duplicate_sequences()
        assert not msa_without_duplicate_sequences.contains_duplicate_sequences()

    def test_get_msa_file_format(self, dna_phylip_msa, dna_fasta_msa):
        ff = dna_phylip_msa._get_file_format()
        assert ff == FileFormat.PHYLIP

        ff = dna_fasta_msa._get_file_format()
        assert ff == FileFormat.FASTA

    def test_get_msa_file_format_raises_value_error(self,
                                                    raxmlng_inference_log):
        with pytest.raises(ValueError):
            MSA(raxmlng_inference_log)._get_file_format()

    def test_guess_msa_file_data_type(self):
        cwd = os.getcwd()
        for true_type in DataType:
            base_dir = os.path.join(cwd, "tests", "data", true_type.value)
            for msa_file in os.listdir(base_dir):
                msa_file = os.path.join(base_dir, msa_file)
                msa = MSA(msa_file)
                guessed_type = msa.guess_data_type()
                assert guessed_type == true_type

    def test_get_raxmlng_model(self, all_msa_files_with_model):
        for msa_file, true_model in all_msa_files_with_model:
            msa = MSA(msa_file)
            model = msa.get_raxmlng_model()

            assert model == true_model

    def test_number_of_taxa(self, dna_phylip_msa):
        assert dna_phylip_msa.number_of_taxa() == 68

    def test_number_of_sites(self, dna_phylip_msa):
        assert dna_phylip_msa.number_of_sites() == 766

    def test_column_entropies(self, dna_phylip_msa):
        column_entropies = dna_phylip_msa.column_entropies()

        assert column_entropies[0] == 0.0
        assert column_entropies[5] == pytest.approx(0.322757, abs=0.01)

    def test_entropy(self, dna_phylip_msa):
        entropy = dna_phylip_msa.entropy()

        assert entropy == pytest.approx(0.19863, abs=0.01)

    def test_bollback_multinomial(self, small_msa):
        bollback = small_msa.bollback_multinomial()
        
        assert bollback == pytest.approx(-365.70127, abs=0.001)

    def test_treelikeness_score(self, dna_phylip_msa):
        return
        # print(dna_phylip_msa.number_of_taxa())
        # print(dna_phylip_msa.treelikeness_score())

    def test_save_reduced_alignment_without_replace_does_not_change_msa(self, msa_with_duplicate_sequences):
        pre_taxa = msa_with_duplicate_sequences.number_of_taxa()
        pre_sites = msa_with_duplicate_sequences.number_of_sites()

        reduced_msa = os.path.join(os.getcwd(), "tests", "data", "DNA", "0.phy.pythia.reduced")

        msa_with_duplicate_sequences.save_reduced_alignment(reduced_msa, replace_original=False)

        post_taxa = msa_with_duplicate_sequences.number_of_taxa()
        post_sites = msa_with_duplicate_sequences.number_of_sites()

        assert pre_taxa == post_taxa
        assert pre_sites == post_sites

    def test_save_reduced_alignment_with_replace_changes_msa(self, msa_with_duplicate_sequences):
        pre_taxa = msa_with_duplicate_sequences.number_of_taxa()
        pre_sites = msa_with_duplicate_sequences.number_of_sites()

        reduced_msa = os.path.join(os.getcwd(), "tests", "data", "DNA", "0.phy.pythia.reduced")
        if os.path.exists(reduced_msa):
            os.remove(reduced_msa)

        msa_with_duplicate_sequences.save_reduced_alignment(reduced_msa, replace_original=True)

        post_taxa = msa_with_duplicate_sequences.number_of_taxa()
        post_sites = msa_with_duplicate_sequences.number_of_sites()

        assert pre_sites == post_sites # number of sites still needs to be the same
        assert pre_taxa > post_taxa

    def test_save_reduced_alignemnt_without_duplicate_raises_pypythia_exception(self, msa_without_duplicate_sequences):
        with pytest.raises(PyPythiaException):
            msa_without_duplicate_sequences.save_reduced_alignment("")
