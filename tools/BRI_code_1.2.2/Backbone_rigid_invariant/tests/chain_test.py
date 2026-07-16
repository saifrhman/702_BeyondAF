# -*- coding = utf-8 -*-
# @File: chain_test.PY

from pathlib import Path

import pytest
import pandas as pd

from bri import Chain, MiniChain

test_data_directory = Path("tests/test_data/")


def local_test_file_str(filename):
    return str(test_data_directory.joinpath(filename))


@pytest.mark.parametrize(
    "pdb_id, entity_id, model_id, chain_id, start_index, chain_length, extra_keys, data",
    [
        ("1HHO", 1, 1, "A", 1, 141, None, None),
        (local_test_file_str("1hho.cif"), 1, 1, "A", 1, 141, None, None),
    ],
)
def test_Chain_functions(
    pdb_id, entity_id, model_id, chain_id, start_index, chain_length, extra_keys, data
):
    chain = Chain(
        pdb_id,
        entity_id,
        model_id,
        chain_id,
        start_index,
        chain_length,
        extra_keys,
        data,
    )

    property_test(
        chain,
        pdb_id,
        entity_id,
        model_id,
        chain_id,
        start_index,
        chain_length,
        extra_keys,
        data,
    )

    invairant_perturb_test(chain)


@pytest.mark.parametrize(
    "pdb_id, model_id, chain_id, start_index, chain_length, extra_keys, data",
    [
        (local_test_file_str("1ilz_extract.cif"), 1, '.', 1, 257, None, None),
    ],
)
def test_MiniChain_functions(
    pdb_id, model_id, chain_id, start_index, chain_length, extra_keys, data
):
    chain = MiniChain(
        pdb_id,
        model_id,
        chain_id,
        start_index,
        chain_length,
        extra_keys,
        data,
    )

    property_test(
        chain,
        pdb_id,
        None,
        model_id,
        chain_id,
        start_index,
        chain_length,
        extra_keys,
        data,
    )

    invairant_perturb_test(chain)


def property_test(
    chain: Chain,
    pdb_id,
    entity_id,
    model_id,
    chain_id,
    start_index,
    chain_length,
    extra_keys,
    data,
):
    assert len(chain) == chain_length
    print(chain.get_chain_invariant_BTP())


def invairant_perturb_test(chain: Chain):
    init_coordinate = chain.coordinates[:]
    init_invariant = chain.invariant[:]

    assert chain.perturb_radius == 0
    assert chain.invariant.empty is False

    chain.perturb_radius = 0.5
    assert chain.perturb_radius == 0.5
    assert init_coordinate.equals(chain.coordinates) is False
    assert init_invariant.equals(chain.invariant) is False

    chain.perturb_radius = 0
    assert chain.perturb_radius == 0
    assert init_coordinate.equals(chain.coordinates) is True
    assert init_invariant.equals(chain.invariant) is True

