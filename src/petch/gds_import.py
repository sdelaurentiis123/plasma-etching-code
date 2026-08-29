"""Minimal, fail-closed GDSII geometry import for feature-mask layouts.

The feature engine needs mask polygons and periodic placement, not a complete
layout-editor implementation.  This module therefore reads the GDSII records
that define polygonal masks (``BOUNDARY``), cell references (``SREF``), and
array references (``AREF``), while retaining the physical database-unit scale.
Unsupported geometric element types are recorded and make a predictive ingest
fail closed unless the caller explicitly accepts them.

GDSII has no material or process semantics.  A parsed polygon is geometry only:
the caller must still identify whether layer polarity means mask solid or mask
opening and must supply the layer stack and surface mechanism separately.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import struct

import numpy as np

__all__ = [
    "GdsArrayReference",
    "GdsCell",
    "GdsLibrary",
    "GdsPolygon",
    "GdsReference",
    "read_gds",
]


def _readonly(value, dtype=float):
    result = np.asarray(value, dtype=dtype).copy()
    result.setflags(write=False)
    return result


def _gds_real8(payload: bytes) -> float:
    """Decode one IBM-base-16 real used by GDSII."""
    if len(payload) != 8:
        raise ValueError("a GDSII real must contain exactly eight bytes")
    if payload == b"\0" * 8:
        return 0.0
    sign = -1.0 if payload[0] & 0x80 else 1.0
    exponent = (payload[0] & 0x7F) - 64
    mantissa = int.from_bytes(payload[1:], "big") / float(1 << 56)
    return sign * mantissa * (16.0 ** exponent)


def _ascii(payload: bytes) -> str:
    return payload.rstrip(b"\0").decode("ascii", errors="strict")


def _int16(payload: bytes) -> tuple[int, ...]:
    if len(payload) % 2:
        raise ValueError("odd-length GDSII int16 payload")
    return struct.unpack(">" + "h" * (len(payload) // 2), payload)


def _xy(payload: bytes) -> np.ndarray:
    if len(payload) % 8:
        raise ValueError("GDSII XY payload is not a sequence of point pairs")
    values = struct.unpack(">" + "i" * (len(payload) // 4), payload)
    return np.asarray(values, dtype=np.int64).reshape(-1, 2)


@dataclass(frozen=True)
class GdsPolygon:
    layer: int
    datatype: int
    vertices_db: np.ndarray

    def __post_init__(self):
        vertices = _readonly(self.vertices_db, dtype=np.int64)
        if vertices.ndim != 2 or vertices.shape[1] != 2 or len(vertices) < 3:
            raise ValueError("a GDSII boundary needs at least three vertices")
        if len(np.unique(vertices, axis=0)) < 3:
            raise ValueError("degenerate GDSII boundary")
        object.__setattr__(self, "vertices_db", vertices)

    @property
    def bounds_db(self) -> tuple[int, int, int, int]:
        lower = self.vertices_db.min(axis=0)
        upper = self.vertices_db.max(axis=0)
        return int(lower[0]), int(lower[1]), int(upper[0]), int(upper[1])


@dataclass(frozen=True)
class GdsReference:
    cell_name: str
    origin_db: np.ndarray
    reflected_x: bool = False
    magnification: float = 1.0
    angle_deg: float = 0.0

    def __post_init__(self):
        origin = _readonly(self.origin_db, dtype=np.int64)
        if origin.shape != (2,) or not self.cell_name:
            raise ValueError("invalid GDSII cell reference")
        if not np.isfinite(self.magnification) or self.magnification <= 0.0:
            raise ValueError("invalid GDSII reference magnification")
        if not np.isfinite(self.angle_deg):
            raise ValueError("invalid GDSII reference angle")
        object.__setattr__(self, "origin_db", origin)

    def transform(self, points_db: np.ndarray) -> np.ndarray:
        points = np.asarray(points_db, dtype=float).copy()
        if self.reflected_x:
            points[:, 1] *= -1.0
        points *= float(self.magnification)
        theta = np.deg2rad(float(self.angle_deg))
        rotation = np.asarray([
            [np.cos(theta), -np.sin(theta)],
            [np.sin(theta), np.cos(theta)],
        ])
        return points @ rotation.T + self.origin_db


@dataclass(frozen=True)
class GdsArrayReference:
    cell_name: str
    columns: int
    rows: int
    origin_db: np.ndarray
    column_vector_db: np.ndarray
    row_vector_db: np.ndarray
    reflected_x: bool = False
    magnification: float = 1.0
    angle_deg: float = 0.0

    def __post_init__(self):
        for name in ("origin_db", "column_vector_db", "row_vector_db"):
            value = _readonly(getattr(self, name), dtype=float)
            if value.shape != (2,) or np.any(~np.isfinite(value)):
                raise ValueError("invalid GDSII array vector")
            object.__setattr__(self, name, value)
        if self.columns <= 0 or self.rows <= 0 or not self.cell_name:
            raise ValueError("invalid GDSII array dimensions")
        if not np.isfinite(self.magnification) or self.magnification <= 0.0:
            raise ValueError("invalid GDSII array magnification")
        if not np.isfinite(self.angle_deg):
            raise ValueError("invalid GDSII array angle")

    @property
    def instance_count(self) -> int:
        return int(self.columns) * int(self.rows)

    def instance_origins_db(self) -> np.ndarray:
        column = np.arange(self.columns, dtype=float)[:, None]
        row = np.arange(self.rows, dtype=float)[:, None]
        return (
            self.origin_db[None, None, :]
            + column[:, None, :] * self.column_vector_db[None, None, :]
            + row[None, :, :] * self.row_vector_db[None, None, :]
        ).reshape(-1, 2)

    def as_reference(self, origin_db) -> GdsReference:
        return GdsReference(
            self.cell_name,
            np.asarray(origin_db, dtype=float),
            reflected_x=self.reflected_x,
            magnification=self.magnification,
            angle_deg=self.angle_deg,
        )


@dataclass(frozen=True)
class GdsCell:
    name: str
    polygons: tuple[GdsPolygon, ...]
    references: tuple[GdsReference, ...]
    arrays: tuple[GdsArrayReference, ...]
    unsupported_element_types: tuple[int, ...] = ()

    @property
    def is_geometry_complete(self) -> bool:
        return not self.unsupported_element_types


@dataclass(frozen=True)
class GdsLibrary:
    name: str
    user_unit_per_database_unit: float
    database_unit_m: float
    cells: dict[str, GdsCell]

    def __post_init__(self):
        if (
            not self.name
            or not np.isfinite(self.user_unit_per_database_unit)
            or self.user_unit_per_database_unit <= 0.0
            or not np.isfinite(self.database_unit_m)
            or self.database_unit_m <= 0.0
            or not self.cells
        ):
            raise ValueError("invalid GDSII library metadata")
        missing = sorted({
            item.cell_name
            for cell in self.cells.values()
            for item in (*cell.references, *cell.arrays)
            if item.cell_name not in self.cells
        })
        if missing:
            raise ValueError(f"GDSII references missing cells: {missing}")

    @property
    def top_cells(self) -> tuple[str, ...]:
        referenced = {
            item.cell_name
            for cell in self.cells.values()
            for item in (*cell.references, *cell.arrays)
        }
        return tuple(sorted(set(self.cells) - referenced))

    @property
    def geometry_complete(self) -> bool:
        return all(cell.is_geometry_complete for cell in self.cells.values())

    def expanded_reference_counts(self, top_cell: str) -> dict[str, int]:
        """Count leaf-cell instances below ``top_cell`` without flattening."""
        if top_cell not in self.cells:
            raise KeyError(top_cell)
        output: dict[str, int] = {}

        def visit(name: str, multiplier: int, ancestry: tuple[str, ...]):
            if name in ancestry:
                raise ValueError(f"cyclic GDSII hierarchy: {ancestry + (name,)}")
            cell = self.cells[name]
            if not cell.references and not cell.arrays:
                output[name] = output.get(name, 0) + int(multiplier)
                return
            for reference in cell.references:
                visit(reference.cell_name, multiplier, ancestry + (name,))
            for array in cell.arrays:
                visit(
                    array.cell_name,
                    multiplier * array.instance_count,
                    ancestry + (name,),
                )

        visit(top_cell, 1, ())
        return dict(sorted(output.items()))


def _records(payload: bytes):
    offset = 0
    while offset < len(payload):
        if offset + 4 > len(payload):
            raise ValueError("truncated GDSII record header")
        length, record_type, data_type = struct.unpack(">HBB", payload[offset:offset + 4])
        if length < 4 or length % 2 or offset + length > len(payload):
            raise ValueError(f"invalid GDSII record length {length} at byte {offset}")
        yield record_type, data_type, payload[offset + 4:offset + length]
        offset += length


def read_gds(source) -> GdsLibrary:
    """Read polygonal GDSII geometry from a path or raw bytes."""
    payload = (
        bytes(source)
        if isinstance(source, (bytes, bytearray))
        else Path(source).read_bytes()
    )
    library_name = None
    unit_pair = None
    cells: dict[str, GdsCell] = {}
    current_cell = None
    current_name = None
    element = None

    for record_type, _data_type, data in _records(payload):
        if record_type == 0x02:  # LIBNAME
            library_name = _ascii(data)
        elif record_type == 0x03:  # UNITS
            if len(data) != 16:
                raise ValueError("invalid GDSII UNITS record")
            unit_pair = (_gds_real8(data[:8]), _gds_real8(data[8:]))
        elif record_type == 0x05:  # BGNSTR
            if current_cell is not None:
                raise ValueError("nested GDSII structure")
            current_cell = {
                "polygons": [], "references": [], "arrays": [],
                "unsupported": [],
            }
            current_name = None
        elif record_type == 0x06:  # STRNAME
            if current_cell is None:
                raise ValueError("GDSII STRNAME outside a structure")
            current_name = _ascii(data)
        elif record_type in (0x08, 0x0A, 0x0B):  # BOUNDARY/SREF/AREF
            if current_cell is None or element is not None:
                raise ValueError("misnested GDSII element")
            element = {"type": record_type}
        elif record_type in (0x09, 0x0C, 0x15, 0x2D):  # PATH/TEXT/NODE/BOX
            if current_cell is None or element is not None:
                raise ValueError("misnested unsupported GDSII element")
            element = {"type": record_type, "unsupported": True}
        elif record_type == 0x0D and element is not None:  # LAYER
            element["layer"] = _int16(data)[0]
        elif record_type == 0x0E and element is not None:  # DATATYPE
            element["datatype"] = _int16(data)[0]
        elif record_type == 0x12 and element is not None:  # SNAME
            element["cell_name"] = _ascii(data)
        elif record_type == 0x13 and element is not None:  # COLROW
            columns, rows = _int16(data)
            element["columns"] = int(columns)
            element["rows"] = int(rows)
        elif record_type == 0x1A and element is not None:  # STRANS
            flags = int.from_bytes(data, "big")
            element["reflected_x"] = bool(flags & 0x8000)
        elif record_type == 0x1B and element is not None:  # MAG
            element["magnification"] = _gds_real8(data)
        elif record_type == 0x1C and element is not None:  # ANGLE
            element["angle_deg"] = _gds_real8(data)
        elif record_type == 0x10 and element is not None:  # XY
            element["xy"] = _xy(data)
        elif record_type == 0x11:  # ENDEL
            if current_cell is None or element is None:
                raise ValueError("GDSII ENDEL without an element")
            kind = element["type"]
            if element.get("unsupported"):
                current_cell["unsupported"].append(kind)
            elif kind == 0x08:
                points = np.asarray(element.get("xy", ()), dtype=np.int64)
                if len(points) >= 2 and np.array_equal(points[0], points[-1]):
                    points = points[:-1]
                current_cell["polygons"].append(GdsPolygon(
                    int(element.get("layer", -1)),
                    int(element.get("datatype", -1)),
                    points,
                ))
            elif kind == 0x0A:
                points = np.asarray(element.get("xy", ()), dtype=np.int64)
                if points.shape != (1, 2):
                    raise ValueError("a GDSII SREF must have exactly one origin")
                current_cell["references"].append(GdsReference(
                    element.get("cell_name", ""),
                    points[0],
                    reflected_x=bool(element.get("reflected_x", False)),
                    magnification=float(element.get("magnification", 1.0)),
                    angle_deg=float(element.get("angle_deg", 0.0)),
                ))
            else:
                points = np.asarray(element.get("xy", ()), dtype=np.int64)
                columns = int(element.get("columns", 0))
                rows = int(element.get("rows", 0))
                if points.shape != (3, 2) or columns <= 0 or rows <= 0:
                    raise ValueError("invalid GDSII AREF geometry")
                # GDSII stores a point one full array span away, so the
                # per-instance displacement divides by N, not N-1.  This is
                # especially important for a one-column/one-row AREF.
                current_cell["arrays"].append(GdsArrayReference(
                    element.get("cell_name", ""), columns, rows, points[0],
                    (points[1] - points[0]) / float(columns),
                    (points[2] - points[0]) / float(rows),
                    reflected_x=bool(element.get("reflected_x", False)),
                    magnification=float(element.get("magnification", 1.0)),
                    angle_deg=float(element.get("angle_deg", 0.0)),
                ))
            element = None
        elif record_type == 0x07:  # ENDSTR
            if current_cell is None or element is not None or not current_name:
                raise ValueError("invalid GDSII structure termination")
            if current_name in cells:
                raise ValueError(f"duplicate GDSII cell {current_name!r}")
            cells[current_name] = GdsCell(
                current_name,
                tuple(current_cell["polygons"]),
                tuple(current_cell["references"]),
                tuple(current_cell["arrays"]),
                tuple(sorted(set(current_cell["unsupported"]))),
            )
            current_cell = None
            current_name = None

    if current_cell is not None or element is not None:
        raise ValueError("unterminated GDSII structure or element")
    if library_name is None or unit_pair is None:
        raise ValueError("GDSII library is missing LIBNAME or UNITS")
    return GdsLibrary(library_name, unit_pair[0], unit_pair[1], cells)
