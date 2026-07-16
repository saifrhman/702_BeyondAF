# -*- coding = utf-8 -*-
# @File: structure_test.PY

from pathlib import Path

import pytest
import pandas as pd

from bri.base.base_util import StructureBase
from bri import Chain, Entry

test_data_directory = Path("tests/test_data/")
a_1hho = pd.read_csv(test_data_directory.joinpath("1HHO-1-1-A-1-141.csv"))


@pytest.mark.parametrize(
    "def_class, init_params",
    [
        (StructureBase, ["1HHO"]),
        (StructureBase, ["1HHO", {"method": "refine.pdbx_refine_id"}]),
        (Chain, ["1HHO", 1, 1, "A", 1, 141]),
        (
            Chain,
            ["1HHO", 1, 1, "A", 1, 141, None, a_1hho],
        ),
        (Chain, ["1HHO", 1, 1, "A", 1, 141, {"method": "refine.pdbx_refine_id"}]),
        (Entry, ["1HHO"]),
    ],
)
def test_base_init(def_class, init_params):
    instance = def_class(*init_params)
    assert instance.pdb_id == Path(init_params[0]).name.upper()
    assert instance.coordinates is not None
    assert list(instance.coordinates.columns).sort() == list(instance.extract_cols.keys()).sort()


def test_StructureBase_functions():
    extra_attr_pass = StructureBase("1HHO", {"method": "refine.pdbx_refine_id"})
    assert extra_attr_pass.method == "X-RAY DIFFRACTION"

    extra_attr_fail = StructureBase("7FFO", {"method": "refine.pdbx_refine_id"})
    assert extra_attr_fail is not None
    assert extra_attr_fail.method is None

    extra_attr_pass = StructureBase(
        "7FFO",
        {
            "exptl_method": "exptl.method",
            "citation": "citation",
            "database": "database_2",
        },
    )
    assert extra_attr_pass.exptl_method == "ELECTRON MICROSCOPY"
    assert isinstance(extra_attr_pass.citation, list)
    print(extra_attr_pass.citation)
    assert isinstance(extra_attr_pass.database, list)
    print(extra_attr_pass.database)
