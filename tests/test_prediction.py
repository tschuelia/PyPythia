import pathlib
import re
import tempfile
import warnings

import pandas as pd
import pytest

from pypythia import __version__
from pypythia.custom_errors import PyPythiaException
from pypythia.msa import parse_msa
from pypythia.prediction import (
    _handle_duplicates,
    _handle_full_gap_sequences,
    collect_features,
    predict_difficulty,
)


def test_handle_duplicates(msa_test_data_row):
    msa = parse_msa(msa_test_data_row.msa_file)
    reduced_msa = _handle_duplicates(msa, deduplicate=True, log_info=False)

    if msa_test_data_row.contains_duplicates:
        assert reduced_msa != msa
        assert reduced_msa.n_taxa < msa.n_taxa
    else:
        assert reduced_msa == msa


def test_handle_duplicates_dont_deduplicate(msa_test_data_row):
    msa = parse_msa(msa_test_data_row.msa_file)
    reduced_msa = _handle_duplicates(msa, deduplicate=False, log_info=False)

    assert reduced_msa == msa


def test_handle_full_gap_sequences(msa_test_data_row):
    msa = parse_msa(msa_test_data_row.msa_file)
    reduced_msa = _handle_full_gap_sequences(msa, remove_full_gaps=True, log_info=False)

    if msa_test_data_row.contains_full_gap_sequences:
        assert reduced_msa != msa
        assert reduced_msa.n_taxa < msa.n_taxa
    else:
        assert reduced_msa == msa


def test_handle_full_gap_sequences_dont_remove_full_gaps(msa_test_data_row):
    msa = parse_msa(msa_test_data_row.msa_file)
    reduced_msa = _handle_full_gap_sequences(
        msa, remove_full_gaps=False, log_info=False
    )

    assert reduced_msa == msa


def test_collect_features(msa_test_data_row, raxmlng):
    msa = parse_msa(msa_test_data_row.msa_file)
    features = collect_features(
        msa=msa, msa_file=msa_test_data_row.msa_file, raxmlng=raxmlng
    )
    assert features.shape[0] == 1

    pd.testing.assert_series_equal(
        features.loc[0],
        msa_test_data_row[features.columns],
        check_dtype=False,
        check_names=False,
    )


def test_collect_features_stores_trees(phylip_msa_file, raxmlng):
    msa = parse_msa(phylip_msa_file)
    with tempfile.NamedTemporaryFile("w") as pars_trees_file:
        pars_trees_file = pathlib.Path(pars_trees_file.name)
        collect_features(
            msa=msa,
            msa_file=phylip_msa_file,
            raxmlng=raxmlng,
            pars_trees_file=pars_trees_file,
        )
        assert pars_trees_file.exists()

        # Should contain 24 parsimony trees
        assert sum(1 for _ in pars_trees_file.open()) == 24


@pytest.mark.parametrize("store_results", [True, False])
@pytest.mark.parametrize("plot_shap", [True, False])
def test_predict_difficulty(
    msa_test_data_row, raxmlng_command, store_results, plot_shap
):
    # Check if the Pythia version is identical, if not the expected difficulty parquet file might be outdated
    # In this case, raise a warning
    if msa_test_data_row.pythia_version != __version__:
        warnings.warn(
            f"The Pythia version in the test data is {msa_test_data_row.pythia_version}, but the current "
            f"Pythia version is {__version__}. The expected difficulty parquet file might be outdated."
        )

    with tempfile.TemporaryDirectory() as tmpdir:
        prefix = pathlib.Path(tmpdir) / "test"

        predicted_difficulty = predict_difficulty(
            msa_file=msa_test_data_row.msa_file,
            raxmlng=raxmlng_command,
            seed=msa_test_data_row.raxmlng_seed,
            deduplicate=False,
            remove_full_gaps=False,
            result_prefix=prefix,
            store_results=store_results,
            plot_shap=plot_shap,
        )

        # 1. Check if the predicted difficulty matches the "ground-truth" in our test data
        assert predicted_difficulty == pytest.approx(
            msa_test_data_row.predicted_difficulty, abs=0.01
        )

        # 2. Check if the results exist if store_results=True, else check if they don't exist
        pars_trees_file = pathlib.Path(f"{prefix}.pythia.trees")
        shap_file = pathlib.Path(f"{prefix}.shap.pdf")
        results_file = pathlib.Path(f"{prefix}.pythia.csv")
        reduced_msa_file = pathlib.Path(f"{prefix}.reduced.phy")

        if store_results:
            assert pars_trees_file.exists()
            assert results_file.exists()
            if plot_shap:
                assert shap_file.exists()
            else:
                assert not shap_file.exists()
        else:
            assert not pars_trees_file.exists()
            assert not shap_file.exists()
            assert not results_file.exists()

        # Since we set deduplicate=False and remove_full_gaps=False, the reduced MSA should be identical to the original
        # and no reduced MSA should be saved in either case
        assert not reduced_msa_file.exists()


def test_predict_difficulty_with_deduplication_and_gap_removal(
    msa_test_data, raxmlng_command
):
    """
    In this test case, we check if the deduplication and gap removal works as expected. For this, we use the test MSAs
    that contain duplicate and/or full-gap sequences and compare the results to the respective data in msa_test_data
    for the reduced MSA.

    If, for instance, the MSA `DNA/test.phy` contains duplicates/full-gap sequences, the respective reduced MSA
    `DNA/test.reduced.phy` should be in `msa_test_data`, and running predict_difficulty with `DNA/test.phy` and
    `deduplicate=True` and `remove_full_gaps=True` should yield the same difficulty as the difficulty for
    `DNA/test.reduced.phy` in `msa_test_data`.
    """

    # Check if the Pythia version is identical, if not the expected difficulty parquet file might be outdated
    # In this case, raise a warning
    pythia_version_test_data = msa_test_data.pythia_version.unique()[0]
    if pythia_version_test_data != __version__:
        warnings.warn(
            f"The Pythia version in the test data is {pythia_version_test_data}, but the current "
            f"Pythia version is {__version__}. The expected difficulty parquet file might be outdated."
        )

    data_with_duplicates_or_full_gaps = msa_test_data.loc[
        lambda x: x.contains_duplicates | x.contains_full_gap_sequences
    ]

    # Sanity check that we actually ran all tests
    expected_n_test_cases = 5
    actual_n_test_cases = 0

    for idx, row in data_with_duplicates_or_full_gaps.iterrows():
        # First, we check if the respective reduced alignment is also in `msa_test_data`
        reduced_msa_file_name = f"{row.data_type}/{row.msa_file.stem}.reduced.phy"
        reduced_row = msa_test_data.loc[
            msa_test_data.msa_file.astype(str).str.endswith(reduced_msa_file_name)
        ]
        if reduced_row.empty:
            continue
        elif reduced_row.shape[0] > 1:
            raise ValueError(
                f"Multiple rows found for reduced MSA {reduced_msa_file_name}"
            )

        with tempfile.TemporaryDirectory() as tmpdir:
            prefix = pathlib.Path(tmpdir) / "test"

            predicted_difficulty = predict_difficulty(
                msa_file=row.msa_file,
                raxmlng=raxmlng_command,
                deduplicate=True,
                remove_full_gaps=True,
                result_prefix=prefix,
                store_results=True,
            )

            # 1. We expect a reduced MSA to be saved
            reduced_msa_file = pathlib.Path(f"{prefix}.reduced.phy")
            assert reduced_msa_file.exists()

            # 2. Check if the predicted difficulty matches the "ground-truth" in our test data
            # Note that we check against the difficulty for the reduced MSA in the test data
            assert predicted_difficulty == pytest.approx(
                reduced_row.predicted_difficulty.mean(), abs=0.01
            )

            # 3. Check if the computed features match the features of the reduced MSA in the test data
            results_file = pathlib.Path(f"{prefix}.pythia.csv")
            results = pd.read_csv(results_file)
            results.drop(
                columns=["difficulty", "msa_file"], inplace=True
            )  # difficulty already with tolerance
            pd.testing.assert_frame_equal(
                results.reset_index(drop=True),
                reduced_row[results.columns].reset_index(drop=True),
            )
            actual_n_test_cases += 1

    assert actual_n_test_cases == expected_n_test_cases, "Not all test cases were run."


def test_predict_difficulty_fewer_than_four_taxa(data_dir, raxmlng_command):
    with pytest.raises(
        PyPythiaException, match="The MSA contains less than 4 sequences."
    ):
        with tempfile.TemporaryDirectory() as tmpdir:
            prefix = pathlib.Path(tmpdir) / "test"
            predict_difficulty(
                msa_file=data_dir / "DNA" / "3_taxa_msa.fasta",
                raxmlng=raxmlng_command,
                deduplicate=False,
                remove_full_gaps=False,
                result_prefix=prefix,
                store_results=True,
            )


def test_predict_difficulty_with_deduplication_and_gap_removal_if_reduced_msa_has_fewer_than_four_taxa(
    data_dir, raxmlng_command
):
    # DNA/5.phy contains 10 taxa, but after deduplication and gap removal, only 3 taxa remain
    # in this case, we expect the general error, but also the hint about the reduced MSA
    with pytest.raises(
        PyPythiaException,
        match=re.compile(
            r"During preprocessing.+leading to an MSA with less than 4 sequences",
            re.DOTALL,
        ),
    ):
        with tempfile.TemporaryDirectory() as tmpdir:
            prefix = pathlib.Path(tmpdir) / "test"
            predict_difficulty(
                msa_file=data_dir / "DNA" / "5.phy",
                raxmlng=raxmlng_command,
                deduplicate=True,
                remove_full_gaps=True,
                result_prefix=prefix,
                store_results=True,
            )


def test_predict_difficulty_msa_files_does_not_exist(raxmlng_command):
    msa_file = pathlib.Path("this_does_not_exist.fasta")
    with pytest.raises(
        PyPythiaException, match=f"The given MSA {msa_file} file does not exist."
    ):
        predict_difficulty(
            msa_file=msa_file,
            raxmlng=raxmlng_command,
            deduplicate=True,
            remove_full_gaps=True,
            result_prefix=None,
            store_results=False,
        )


def test_predict_difficulty_raxmlng_executable_none(small_msa_file):
    with pytest.raises(
        PyPythiaException, match="Path to the RAxML-NG executable is required"
    ):
        predict_difficulty(
            msa_file=small_msa_file,
            raxmlng=None,
            deduplicate=False,
            remove_full_gaps=False,
            result_prefix=None,
            store_results=False,
        )


def test_predict_difficulty_raxmlng_init_fails(small_msa_file):
    with pytest.raises(PyPythiaException, match="Initializing RAxML-NG failed"):
        predict_difficulty(
            msa_file=small_msa_file,
            raxmlng=pathlib.Path("/path/to/non/existing/raxml-ng"),
            deduplicate=False,
            remove_full_gaps=False,
            result_prefix=None,
            store_results=False,
        )


def test_predict_difficulty_predictor_init_fails(small_msa_file, raxmlng_command):
    with pytest.raises(
        PyPythiaException, match="Initializing the difficulty predictor failed"
    ):
        with tempfile.NamedTemporaryFile("w") as f:
            f.write("This is not a valid predictor file.")
            predict_difficulty(
                msa_file=small_msa_file,
                raxmlng=raxmlng_command,
                deduplicate=False,
                remove_full_gaps=False,
                result_prefix=None,
                store_results=False,
                model_file=pathlib.Path(f.name),
            )
