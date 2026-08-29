from scripts.audit_zhu_npg80_gds_geometry import build


def test_supplied_zhu_gds_is_checksum_bound_and_physically_decoded():
    audit = build()

    assert audit["source"]["sha256"] == (
        "1378b31c6b206a5b62c0254979c0fce36219dc625f716cbd00f858b53f36832b"
    )
    assert audit["gdsii"]["library_name"] == "DOSE_TEST_700_ZEP"
    assert audit["gdsii"]["database_unit_nm"] == 0.5
    assert audit["gdsii"]["unique_array_pitch_nm"] == {"x": 350.0, "y": 350.0}
    assert audit["gdsii"]["top_aref_count"] == 16600
    assert audit["gdsii"]["top_sref_count"] == 85693
    assert len(audit["exact_mask_primitives"]["square_cells"]) == 23
    assert audit["simulation_board"]["square_widths_nm"] == [
        105.0, 110.0, 115.0, 120.0, 125.0, 130.0, 135.0, 140.0,
        145.0, 150.0, 155.0, 160.0, 165.0, 170.0, 185.0, 190.0,
        195.0, 225.0, 230.0, 235.0, 240.0, 245.0, 250.0,
    ]
    assert audit["simulation_board"]["special_cells"] == [
        "CROSS_250x105", "INVHOLE_105", "INVHOLE_250", "RECT_250x105",
    ]
    assert audit["interpretation_limits"]["mask_polarity_confirmed_by_operator"] is False
