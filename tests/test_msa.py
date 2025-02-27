import pathlib
import tempfile

import numpy as np
import pytest

from pypythia.custom_errors import PyPythiaException
from pypythia.custom_types import DataType, FileFormat
from pypythia.msa import (
    MSA,
    _get_file_format,
    _guess_dtype,
    deduplicate_sequences,
    parse,
    remove_full_gap_sequences,
)


def test_parse(msa_test_data):
    for idx, row in msa_test_data.iterrows():
        msa_file = pathlib.Path(row.msa_file)
        msa = parse(msa_file)
        assert msa.n_taxa == row.n_taxa
        assert msa.n_sites == row.n_sites
        assert msa.data_type.value == row.data_type
        assert msa.name == msa_file.name


def test_parse_large_phylip(phylip_msa_file):
    msa = parse(phylip_msa_file)
    assert msa.n_taxa == 68
    assert msa.n_sites == 766
    assert msa.data_type == DataType.DNA
    assert msa.name == phylip_msa_file.name


def test_parse_small_fasta(small_msa_file):
    msa = parse(small_msa_file)
    assert msa.n_taxa == 10
    assert msa.n_sites == 522
    assert msa.data_type == DataType.DNA
    assert msa.name == small_msa_file.name

    np.testing.assert_array_equal(msa.taxa, [f"TAXON{i}" for i in range(10)])
    np.testing.assert_array_equal(
        np.sort(np.unique(msa.sequences)),
        np.sort(np.array([b"A", b"C", b"G", b"T", b"-"])),
    )


def test_msa_init():
    taxa = np.array(["TAXON1", "TAXON2", "TAXON3"])
    sequences = np.array(
        [
            [b"A", b"C", b"G", b"T"],
            [b"A", b"C", b"G", b"T"],
            [b"A", b"C", b"G", b"T"],
        ]
    )

    msa = MSA(taxa, sequences, DataType.DNA, "test")
    assert msa.n_taxa == 3
    assert msa.n_sites == 4
    assert msa.data_type == DataType.DNA
    assert msa.name == "test"

    np.testing.assert_array_equal(msa.taxa, taxa)
    np.testing.assert_array_equal(msa.sequences, sequences)


def test_msa_init_wrong_taxa():
    taxa = np.array(["TAXON1", "TAXON2"])
    sequences = np.array(
        [
            [b"A", b"C", b"G", b"T"],
            [b"A", b"C", b"G", b"T"],
            [b"A", b"C", b"G", b"T"],
        ]
    )

    with pytest.raises(
        PyPythiaException, match="Number of taxa and sequences do not match"
    ):
        MSA(taxa, sequences, DataType.DNA, "test")


def test_msa_str_and_repr():
    taxa = np.array(["TAXON1", "TAXON2", "TAXON3"])
    sequences = np.array(
        [
            [b"A", b"C", b"G", b"T"],
            [b"A", b"C", b"G", b"T"],
            [b"A", b"C", b"G", b"T"],
        ]
    )

    msa = MSA(taxa, sequences, DataType.DNA, "test")
    expected_str = "MSA(name=test, n_taxa=3, n_sites=4, data_type=DNA)"
    assert str(msa) == expected_str
    assert repr(msa) == expected_str


def test_contains_duplicate_sequences(msa_test_data):
    for idx, row in msa_test_data.iterrows():
        msa_file = pathlib.Path(row.msa_file)
        msa = parse(msa_file)
        assert msa.contains_duplicate_sequences() == row.contains_duplicates


def test_contains_full_gap_sequences(msa_test_data):
    for idx, row in msa_test_data.iterrows():
        msa_file = pathlib.Path(row.msa_file)
        msa = parse(msa_file)
        assert msa.contains_full_gap_sequences() == row.contains_full_gap_sequences


def test_get_file_format(msa_test_data):
    for idx, row in msa_test_data.iterrows():
        msa_file = pathlib.Path(row.msa_file)
        assert _get_file_format(msa_file).value == row.file_format


def test_get_msa_file_format_raises_value_error(raxmlng_inference_log):
    with pytest.raises(
        PyPythiaException,
        match=f"The file type of {raxmlng_inference_log} could not be determined.",
    ):
        _get_file_format(raxmlng_inference_log)


def test_guess_dtype(msa_test_data):
    for idx, row in msa_test_data.iterrows():
        msa_file = pathlib.Path(row.msa_file)
        msa = parse(msa_file)
        assert _guess_dtype(msa.sequences).value == row.data_type


def test_guess_dtype_fails():
    sequences = np.array(
        [
            [b"A", b"C", b"G", b"T"],
            [b"A", b"C", b"G", b"T"],
            [
                b"J",
                b"C",
                b"G",
                b"T",
            ],  # contains the character J which is neither DNA nor AA
        ]
    )

    with pytest.raises(
        PyPythiaException, match="Data type for character set could not be inferred"
    ):
        _guess_dtype(sequences)


def test_get_raxmlng_model(msa_test_data):
    for idx, row in msa_test_data.iterrows():
        msa_file = pathlib.Path(row.msa_file)
        msa = parse(msa_file)
        assert msa.get_raxmlng_model() == row.raxmlng_model


def test_get_raxmlng_model_fails_for_invalid_dtype():
    with pytest.raises(PyPythiaException, match="Unsupported data type:"):
        MSA(
            np.array(["TAXON1", "TAXON2", "TAXON3"]),
            np.array(
                [
                    [b"A", b"C", b"G", b"T"],
                    [b"A", b"C", b"G", b"T"],
                    [b"A", b"C", b"G", b"T"],
                ]
            ),
            "INVALID_DTYPE",
            "test",
        ).get_raxmlng_model()


def test_write(phylip_msa_file):
    msa = parse(phylip_msa_file)
    with tempfile.NamedTemporaryFile() as tmpfile:
        tmpfile = pathlib.Path(tmpfile.name)
        msa.write(tmpfile, file_format=FileFormat.PHYLIP)

        # File format is correct
        assert _get_file_format(tmpfile) == FileFormat.PHYLIP

        # Number of taxa and sites is identical
        msa_reread = parse(tmpfile)
        assert msa_reread.n_taxa == msa.n_taxa
        assert msa_reread.n_sites == msa.n_sites

        # Names of taxa are identical
        np.testing.assert_array_equal(msa_reread.taxa, msa.taxa)

        # Sequences are identical
        np.testing.assert_array_equal(msa_reread.sequences, msa.sequences)

        # Data type is identical
        assert msa_reread.data_type == msa.data_type


class TestMSAFeatures:
    def test_n_taxa(self, msa_test_data):
        for idx, row in msa_test_data.iterrows():
            msa_file = pathlib.Path(row.msa_file)
            msa = parse(msa_file)
            assert msa.n_taxa == row.n_taxa

    def test_n_sites(self, msa_test_data):
        for idx, row in msa_test_data.iterrows():
            msa_file = pathlib.Path(row.msa_file)
            msa = parse(msa_file)
            assert msa.n_sites == row.n_sites

    def test_n_patterns(self, msa_test_data):
        for idx, row in msa_test_data.iterrows():
            msa_file = pathlib.Path(row.msa_file)
            msa = parse(msa_file)
            assert msa.n_patterns == row.n_patterns

    def test_percentage_gaps(self, msa_test_data):
        for idx, row in msa_test_data.iterrows():
            msa_file = pathlib.Path(row.msa_file)
            msa = parse(msa_file)
            assert msa.proportion_gaps == row.proportion_gaps

    def test_percentage_invariant(self, msa_test_data):
        for idx, row in msa_test_data.iterrows():
            msa_file = pathlib.Path(row.msa_file)
            msa = parse(msa_file)
            assert msa.proportion_invariant == row.proportion_invariant

    def test_entropy(self, msa_test_data):
        for idx, row in msa_test_data.iterrows():
            msa_file = pathlib.Path(row.msa_file)
            msa = parse(msa_file)
            assert msa.entropy() == row.entropy

    def test_pattern_entropy(self, msa_test_data):
        for idx, row in msa_test_data.iterrows():
            msa_file = pathlib.Path(row.msa_file)
            msa = parse(msa_file)
            assert msa.pattern_entropy() == row.pattern_entropy

    def test_bollback_multinomial(self, msa_test_data):
        for idx, row in msa_test_data.iterrows():
            msa_file = pathlib.Path(row.msa_file)
            msa = parse(msa_file)
            assert msa.bollback_multinomial() == row.bollback_multinomial


def test_remove_full_gap_sequences(msa_test_data):
    for idx, row in msa_test_data.iterrows():
        msa_file = pathlib.Path(row.msa_file)
        msa = parse(msa_file)

        if row.contains_full_gap_sequences:
            # If the MSA contains full-gap sequences: expect these sequences to be removed
            msa_no_full_gaps = remove_full_gap_sequences(msa)
            assert not msa_no_full_gaps.contains_full_gap_sequences()
            assert msa_no_full_gaps.n_taxa < msa.n_taxa
            # Number of sites should not be affected
            assert msa_no_full_gaps.n_sites == msa.n_sites
        else:
            # Otherwise, expect a PyPythiaException
            with pytest.raises(
                PyPythiaException, match="No full-gap sequences found in MSA."
            ):
                remove_full_gap_sequences(msa)


def test_deduplicate_sequences(msa_test_data):
    for idx, row in msa_test_data.iterrows():
        msa_file = pathlib.Path(row.msa_file)
        msa = parse(msa_file)

        if row.contains_duplicates:
            # If the MSA contains duplicate sequences: expect these sequences to be removed
            msa_no_duplicates = deduplicate_sequences(msa)
            assert not msa_no_duplicates.contains_duplicate_sequences()
            assert msa_no_duplicates.n_taxa < msa.n_taxa
            # Number of sites should not be affected
            assert msa_no_duplicates.n_sites == msa.n_sites
        else:
            # Otherwise, expect a PyPythiaException
            with pytest.raises(
                PyPythiaException, match="No duplicate sequences found in MSA."
            ):
                deduplicate_sequences(msa)
