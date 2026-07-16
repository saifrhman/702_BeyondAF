# -*- coding = utf-8 -*-
# @File: clean_test.PY

from pathlib import Path

import pytest
import pandas as pd

from bri import Entry
from bri.filter import entry_integrated_cleaning, minientry_integrated_cleaning

test_data_directory = Path("tests/test_data/")


def local_test_file_str(filename):
    return str(test_data_directory.joinpath(filename))


@pytest.mark.parametrize(
    "pdb_path",
    [
        ("1HHO"),
        ("2HHB"),
        (test_data_directory / "AF-G5EB01-F1-model_v4.cif"),
    ],
)
def test_entry_integrated_clean(pdb_path):
    clean_set, dirty_set = entry_integrated_cleaning(pdb_path)
    print(dirty_set[["chain_id", "type"]])

    assert dirty_set.empty is True


@pytest.mark.parametrize(
    "pdb_path",
    [
        (local_test_file_str("1ilz_extract.cif")),
        (local_test_file_str("AF-G5EB01-F1-model_v4.cif")),
        (local_test_file_str("AF-A0A0G3QJE1-F1-model_v4_sel80.0_extract.cif")),
        ("AF_AFA0A023PZF5F1"),
        ("MA_MAT3VR3952"),
    ],
)
def test_minientry_integrated_clean(pdb_path):
    clean_set, dirty_set = minientry_integrated_cleaning(pdb_path)

    assert dirty_set.empty is True
