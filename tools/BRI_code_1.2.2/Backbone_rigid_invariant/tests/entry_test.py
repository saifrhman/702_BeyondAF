# -*- coding = utf-8 -*-
# @File: entry_test.PY

from pathlib import Path

import pytest
import pandas as pd

from bri import Entry, MiniEntry

test_data_directory = Path("tests/test_data/")
a_1hho = pd.read_csv(test_data_directory.joinpath("1HHO-1-1-A-1-141.csv"))


def local_test_file_str(filename: str):
    return str(test_data_directory.joinpath(filename))


@pytest.mark.parametrize(
    "path, extra_keys, data",
    [
        ("1HHO", None, None),
        ("2HHB", None, None),
        ("7T3M", None, None),
    ],
)
def test_Entry_functions(path, extra_keys, data):
    entry = Entry(path, extra_keys, data)

    # basic info
    property_test(entry, extra_keys, data)

    # invariant
    invairant_perturb_test(entry)

    # TODO value check


def property_test(
    entry: Entry,
    extra_keys,
    data,
):
    assert entry.invariant is not None


def invairant_perturb_test(entry: Entry):
    init_coordinate = entry.coordinates[:]
    init_invariant = entry.invariant[:]

    assert entry.perturb_radius == 0
    assert entry.invariant.empty is False

    entry.perturb_radius = 0.5
    assert entry.perturb_radius == 0.5
    assert init_coordinate.equals(entry.coordinates) is False
    assert init_invariant.equals(entry.invariant) is False

    entry.perturb_radius = 0
    assert entry.perturb_radius == 0
    assert init_coordinate.equals(entry.coordinates) is True
    assert init_invariant.equals(entry.invariant) is True


@pytest.mark.parametrize(
    "path, extra_keys, data",
    [
        (local_test_file_str("AF-G5EB01-F1-model_v4.cif"), None, None),
        (
            local_test_file_str("AF-A0A0G3QJE1-F1-model_v4_sel80.0_extract.cif"),
            None,
            None,
        ),
        ("AF_AFA0A023PZF5F1", None, None),
    ],
)
def test_MiniEntry_functions(path, extra_keys, data):
    entry = MiniEntry(path, extra_keys, data)

    # basic info
    property_test(entry, extra_keys, data)

    # invariant
    invairant_perturb_test(entry)
