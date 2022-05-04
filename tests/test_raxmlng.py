from fixtures import *
from tempfile import TemporaryDirectory


def test_infer_parsimony_trees(raxmlng, example_msa_path):
    with TemporaryDirectory() as tmpdir:
        file_path = raxmlng.infer_parsimony_trees(
            msa_file=example_msa_path,
            model="GTR+G",
            prefix=tmpdir,
            n_trees=10
        )

        trees = open(file_path).readlines()
        trees = [t.strip() for t in trees if t]

        assert len(trees) == 10


def test_get_rfdistance_results(raxmlng, multiple_trees_path):
    num_topos, rel_rfdist, abs_rfdist = raxmlng.get_rfdistance_results(multiple_trees_path)
    assert num_topos == 6
    assert rel_rfdist == pytest.approx(0.114, abs=0.01)
    assert abs_rfdist == pytest.approx(22.269, abs=0.1)


def test_get_patterns_gaps_invariant(raxmlng, example_msa_path):
    patterns, gaps, invariant = raxmlng.get_patterns_gaps_invariant(example_msa_path, "GTR+G")

    assert patterns == 241
    assert gaps == pytest.approx(0.0789)
    assert invariant == pytest.approx(0.5979)