"""Gates for deterministic GDSII polygon and hierarchy ingest."""

import struct

import numpy as np
import pytest

from petch.gds_import import read_gds


def _record(record_type, data_type=0, payload=b""):
    length = 4 + len(payload)
    assert length % 2 == 0
    return struct.pack(">HBB", length, record_type, data_type) + payload


def _string(value):
    payload = value.encode("ascii")
    return payload if len(payload) % 2 == 0 else payload + b"\0"


def _real8(value):
    if value == 0.0:
        return b"\0" * 8
    sign = 0x80 if value < 0.0 else 0
    value = abs(float(value))
    exponent = 64
    while value >= 1.0:
        value /= 16.0
        exponent += 1
    while value < 1.0 / 16.0:
        value *= 16.0
        exponent -= 1
    mantissa = int(round(value * (1 << 56)))
    return bytes([sign | exponent]) + mantissa.to_bytes(7, "big")


def _xy(points):
    flat = [coordinate for point in points for coordinate in point]
    return struct.pack(">" + "i" * len(flat), *flat)


def _minimal_array_library():
    timestamp = struct.pack(">12h", *([0] * 12))
    square = [(-250, -250), (250, -250), (250, 250), (-250, 250), (-250, -250)]
    return b"".join((
        _record(0x00, 2, struct.pack(">h", 600)),
        _record(0x01, 2, timestamp),
        _record(0x02, 6, _string("TEST")),
        _record(0x03, 5, _real8(5e-4) + _real8(5e-10)),
        _record(0x05, 2, timestamp),
        _record(0x06, 6, _string("S250")),
        _record(0x08),
        _record(0x0D, 2, struct.pack(">h", 1)),
        _record(0x0E, 2, struct.pack(">h", 0)),
        _record(0x10, 3, _xy(square)),
        _record(0x11),
        _record(0x07),
        _record(0x05, 2, timestamp),
        _record(0x06, 6, _string("main")),
        _record(0x0B),
        _record(0x12, 6, _string("S250")),
        _record(0x13, 2, struct.pack(">2h", 2, 3)),
        _record(0x10, 3, _xy(((0, 0), (1400, 0), (0, 2100)))),
        _record(0x11),
        _record(0x07),
        _record(0x04),
    ))


def test_gds_reader_preserves_physical_units_polygons_and_array_pitch():
    library = read_gds(_minimal_array_library())

    assert library.name == "TEST"
    assert library.database_unit_m == pytest.approx(0.5e-9)
    assert library.top_cells == ("main",)
    assert library.geometry_complete
    square = library.cells["S250"].polygons[0]
    assert square.layer == 1 and square.datatype == 0
    assert square.bounds_db == (-250, -250, 250, 250)

    array = library.cells["main"].arrays[0]
    assert array.instance_count == 6
    assert np.array_equal(array.column_vector_db, [700.0, 0.0])
    assert np.array_equal(array.row_vector_db, [0.0, 700.0])
    origins = array.instance_origins_db()
    assert set(map(tuple, origins)) == {
        (0.0, 0.0), (0.0, 700.0), (0.0, 1400.0),
        (700.0, 0.0), (700.0, 700.0), (700.0, 1400.0),
    }
    assert library.expanded_reference_counts("main") == {"S250": 6}


def test_missing_reference_and_truncated_record_fail_closed():
    payload = _minimal_array_library().replace(b"S250", b"NOPE", 1)
    with pytest.raises(ValueError, match="references missing cells"):
        read_gds(payload)
    with pytest.raises(ValueError, match="truncated|invalid"):
        read_gds(_minimal_array_library()[:-1])


def test_zero_block_padding_after_endlib_is_accepted_but_nonzero_tail_is_not():
    payload = _minimal_array_library()

    padded = read_gds(payload + b"\0" * 1220)
    assert padded.name == "TEST"
    assert padded.cells["S250"].polygons[0].bounds_db == (
        -250, -250, 250, 250)

    with pytest.raises(ValueError, match="nonzero data after GDSII ENDLIB"):
        read_gds(payload + b"\0" * 31 + b"x")
