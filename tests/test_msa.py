from fixtures import *
from msa import MSA


class TestMSAFeatures:
    def test_get_msa_file_format(self, dna_phylip_msa, dna_fasta_msa):
        ff = dna_phylip_msa._get_file_format()
        assert ff == "phylip-relaxed"

        ff = dna_fasta_msa._get_file_format()
        assert ff == "fasta"

    def test_get_msa_file_format_raises_value_error(self,
                                                    raxmlng_inference_log):
        with pytest.raises(ValueError):
            MSA(raxmlng_inference_log)._get_file_format()

    def test_guess_msa_file_data_type(self):
        cwd = os.getcwd()
        for true_type in ["DNA", "AA"]:
            base_dir = os.path.join(cwd, "tests", "data", true_type)
            for msa_file in os.listdir(base_dir):
                msa_file = os.path.join(base_dir, msa_file)
                msa = MSA(msa_file)
                guessed_type = msa.guess_data_type()
                assert guessed_type == true_type

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
