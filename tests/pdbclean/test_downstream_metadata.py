import gzip

from pdbclean.downstream_metadata import (
    parse_entry_metadata_bytes,
)


def _gzip(text: str) -> bytes:
    return gzip.compress(
        text.encode("utf-8")
    )


def test_parse_ordinary_xray_entry_metadata():
    metadata = parse_entry_metadata_bytes(
        _gzip(
            """
data_test
_exptl.method 'X-RAY DIFFRACTION'
_refine.ls_d_res_high 2.20
_pdbx_database_status.recvd_initial_deposition_date 1997-08-14
_struct.title 'ordinary protein structure'
_struct_keywords.text 'TRANSFERASE'
"""
        ),
        pdb_id="10gs",
    )

    assert metadata.pdb_id == "10gs"
    assert metadata.experimental_methods == (
        "X-RAY DIFFRACTION",
    )
    assert metadata.refine_ls_d_res_high == (
        2.20,
    )
    assert (
        metadata.em_3d_reconstruction_resolution
        == ()
    )
    assert metadata.has_deposit_group is False
    assert (
        metadata.deposit_group_mentions_pandda
        is False
    )
    assert metadata.entry_mentions_pandda is False


def test_parse_pandda_deposit_group_metadata():
    metadata = parse_entry_metadata_bytes(
        _gzip(
            """
data_test
_exptl.method 'X-RAY DIFFRACTION'
_refine.ls_d_res_high 1.87
_struct.title
'PanDDA analysis group deposition -- example'
_struct_keywords.text
'PanDDA fragment screening'
_pdbx_deposit_group.group_id G_1002023
_pdbx_deposit_group.group_title
'PanDDA analysis group deposition of models'
_pdbx_deposit_group.group_description
'PanDDA event maps'
_pdbx_deposit_group.group_type 'ground state'
"""
        ),
        pdb_id="5pp5",
    )

    assert metadata.has_deposit_group is True
    assert (
        metadata.deposit_group_mentions_pandda
        is True
    )
    assert metadata.entry_mentions_pandda is True
    assert metadata.deposit_group_ids == (
        "G_1002023",
    )
    assert metadata.deposit_group_types == (
        "ground state",
    )


def test_preserves_multiple_resolution_values_without_filtering():
    metadata = parse_entry_metadata_bytes(
        _gzip(
            """
data_test
loop_
_refine.entry_id
_refine.ls_d_res_high
TEST 2.50
TEST 3.00
loop_
_em_3d_reconstruction.entry_id
_em_3d_reconstruction.resolution
TEST 3.20
TEST 3.40
"""
        ),
        pdb_id="test",
    )

    assert metadata.refine_ls_d_res_high == (
        2.5,
        3.0,
    )
    assert (
        metadata.em_3d_reconstruction_resolution
        == (
            3.2,
            3.4,
        )
    )
