"""STL geometry import: triangle soup -> certified level-set feature geometry.

Two consumers, one importer:

* **general 3-D geometry** — a watertight STL is rasterized to a signed distance
  field on a petch grid and handed to the common feature engine as a
  :class:`~petch.feature_geometry_state_3d.FeatureGeometry3D`;
* **high-aspect-ratio holes** — a body of revolution is detected, its generator
  profile ``r(z)`` extracted with a deviation receipt, and routed to the
  Clausing-validated axisymmetric exchange operator
  (:mod:`petch.axisymmetric_exchange_3d`) instead of paying for a 3-D grid.

Nothing here is fitted or sampled: the distance is the exact point-triangle
distance to the triangle soup, and the inside/outside decision is the
Van Oosterom-Strackee solid-angle winding number, which is exact for a closed
consistently-oriented surface and has no ray-degeneracy failure mode at
vertices, edges, or coplanar facets.

Sign convention
---------------
petch stores ``phi > 0`` inside solid and ``phi < 0`` in gas (the opposite of
the common graphics convention).  ``rasterize_signed_distance`` returns the
petch convention directly; ``solid_region`` selects whether the STL interior is
the solid ("interior", a mask/substrate body) or the void ("exterior", where the
STL is the cavity/hole body and the solid is everything else in the box).

Resolution requirements
-----------------------
The distance field is exact at the grid nodes, but everything downstream sees
only its zero level set through marching cubes, so the STL is resolved to the
grid, not to the file:

* a feature of width ``w`` needs ``w / dx >= 2`` to exist at all and
  ``w / dx >= 6``-``8`` before its width converges to O(dx^2);
* sharp edges and corners are rounded at the ``dx`` scale (a level set has no
  subcell corner representation);
* a wall thinner than ``2 dx`` is not representable and will either vanish or
  weld shut -- ``rasterize_signed_distance`` reports the narrowest resolved
  band it can see (``min_feature_cells``) so the caller can refuse rather than
  silently etch a different structure.

The faceting of the STL itself is a second, independent error: a cylinder
tessellated with ``n`` azimuthal facets has radius error up to
``R (1 - cos(pi / n))``, which is reported by the axisymmetric extractor as the
measured deviation rather than assumed away.
"""

from __future__ import annotations

import struct
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .axisymmetric_exchange_3d import AxisymmetricProfile
from .feature_geometry_state_3d import FeatureGeometry3D

__all__ = [
    "StlMesh",
    "MeshDiagnostics",
    "SdfRasterReport",
    "AxisymmetryReport",
    "read_stl",
    "write_stl",
    "diagnose_mesh",
    "drop_degenerate_faces",
    "revolved_stl_mesh",
    "rasterize_signed_distance",
    "extract_axisymmetric_profile",
    "to_axisymmetric_profile",
    "assign_materials_by_z",
    "build_feature_geometry_from_stl",
]

_WELD_RELATIVE_TOLERANCE = 1e-9


def _readonly(value, dtype=float):
    result = np.asarray(value, dtype=dtype).copy()
    result.setflags(write=False)
    return result


@dataclass(frozen=True)
class StlMesh:
    """Indexed triangle mesh in file units (welded vertices, shared topology).

    ``file_normals`` retains the per-facet normals as stored by the writer; they
    are advisory only.  Every geometric quantity below is recomputed from the
    vertices, because STL writers routinely emit zero or inconsistent normals.
    """

    vertices: np.ndarray
    faces: np.ndarray
    file_normals: np.ndarray | None = None

    def __post_init__(self):
        vertices = _readonly(self.vertices)
        faces = _readonly(self.faces, dtype=int)
        if (vertices.ndim != 2 or vertices.shape[1] != 3 or len(vertices) < 3
                or faces.ndim != 2 or faces.shape[1] != 3 or len(faces) < 1
                or np.any(~np.isfinite(vertices))
                or np.any(faces < 0) or np.any(faces >= len(vertices))):
            raise ValueError("invalid STL triangle mesh")
        normals = (None if self.file_normals is None
                   else _readonly(self.file_normals))
        if normals is not None and normals.shape != faces.shape:
            raise ValueError("facet normal count does not match the faces")
        object.__setattr__(self, "vertices", vertices)
        object.__setattr__(self, "faces", faces)
        object.__setattr__(self, "file_normals", normals)

    @property
    def triangles(self):
        """(F, 3, 3) triangle vertex coordinates."""
        return self.vertices[self.faces]

    @property
    def face_areas(self):
        tri = self.triangles
        cross = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
        return 0.5 * np.linalg.norm(cross, axis=1)

    @property
    def face_normals(self):
        """Unit normals recomputed from the vertex winding (zero if degenerate)."""
        tri = self.triangles
        cross = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
        norm = np.linalg.norm(cross, axis=1)
        safe = np.where(norm > 0.0, norm, 1.0)[:, None]
        return np.where(norm[:, None] > 0.0, cross / safe, 0.0)

    @property
    def bounds(self):
        """((xmin, ymin, zmin), (xmax, ymax, zmax))."""
        return (tuple(self.vertices.min(axis=0)), tuple(self.vertices.max(axis=0)))

    @property
    def signed_volume(self):
        """Divergence-theorem volume; positive when facets wind outward."""
        tri = self.triangles
        return float(np.sum(np.einsum(
            "fj,fj->f", tri[:, 0], np.cross(tri[:, 1], tri[:, 2]))) / 6.0)

    def flipped(self):
        """Same surface with every facet winding reversed."""
        return StlMesh(
            self.vertices, self.faces[:, ::-1],
            None if self.file_normals is None else -self.file_normals)


def _weld(triangle_vertices, relative_tolerance=_WELD_RELATIVE_TOLERANCE):
    """Merge coincident corner coordinates into a shared-vertex topology.

    STL stores every facet independently, so watertightness is only decidable
    after welding.  The quantum is relative to the model extent, not absolute,
    so the same mesh welds identically in nm and in m.
    """
    corners = np.asarray(triangle_vertices, dtype=float).reshape(-1, 3)
    extent = float(np.max(corners.max(axis=0) - corners.min(axis=0)))
    quantum = relative_tolerance * (extent if extent > 0.0 else 1.0)
    keys = np.round(corners / quantum).astype(np.int64)
    _, first, inverse = np.unique(keys, axis=0, return_index=True, return_inverse=True)
    vertices = corners[first]
    faces = inverse.reshape(-1, 3)
    return vertices, faces


def _read_binary_stl(payload):
    count = struct.unpack("<I", payload[80:84])[0]
    expected = 84 + 50 * count
    if len(payload) < expected:
        raise ValueError(
            f"truncated binary STL: {count} facets declared, "
            f"{(len(payload) - 84) // 50} present")
    records = np.frombuffer(
        payload[84:expected], dtype=np.dtype([
            ("normal", "<f4", 3), ("corners", "<f4", (3, 3)),
            ("attribute", "<u2")]), count=count)
    return (np.asarray(records["corners"], dtype=float),
            np.asarray(records["normal"], dtype=float))


def _read_ascii_stl(text):
    normals = []
    corners = []
    current = []
    for line in text.splitlines():
        fields = line.split()
        if not fields:
            continue
        keyword = fields[0].lower()
        if keyword == "facet" and len(fields) >= 5:
            normals.append([float(value) for value in fields[-3:]])
        elif keyword == "vertex" and len(fields) >= 4:
            current.append([float(value) for value in fields[-3:]])
        elif keyword == "endfacet":
            if len(current) != 3:
                raise ValueError(
                    f"ASCII STL facet has {len(current)} vertices, expected 3")
            corners.append(current)
            current = []
    if not corners:
        raise ValueError("ASCII STL contains no facets")
    if len(normals) != len(corners):
        normals = [[0.0, 0.0, 0.0]] * len(corners)
    return np.asarray(corners, dtype=float), np.asarray(normals, dtype=float)


def _looks_binary(payload):
    """Binary iff the declared facet count matches the byte length exactly.

    The usual "starts with solid" sniff is unreliable -- binary writers put
    arbitrary text in the 80-byte header, including "solid".  The length
    relation is decisive.
    """
    if len(payload) < 84:
        return False
    count = struct.unpack("<I", payload[80:84])[0]
    return len(payload) == 84 + 50 * count


def read_stl(source, *, weld_relative_tolerance=_WELD_RELATIVE_TOLERANCE):
    """Read a binary or ASCII STL into a welded :class:`StlMesh`.

    ``source`` is a path or raw ``bytes``.  The format is detected from the
    byte-length relation, not from the leading keyword.
    """
    payload = (bytes(source) if isinstance(source, (bytes, bytearray))
               else Path(source).read_bytes())
    if _looks_binary(payload):
        corners, normals = _read_binary_stl(payload)
    else:
        corners, normals = _read_ascii_stl(payload.decode("utf-8", errors="replace"))
    vertices, faces = _weld(corners, weld_relative_tolerance)
    return StlMesh(vertices, faces, normals)


def write_stl(path, mesh, *, binary=True, name="petch"):
    """Write ``mesh`` as STL (binary by default); returns ``path``."""
    tri = np.asarray(mesh.triangles, dtype=float)
    normals = np.asarray(mesh.face_normals, dtype=float)
    path = Path(path)
    if binary:
        payload = bytearray(struct.pack("<80sI", name.encode()[:80], len(tri)))
        for facet, normal in zip(tri, normals):
            payload += struct.pack("<12fH", *normal, *facet.reshape(-1), 0)
        path.write_bytes(bytes(payload))
        return path
    lines = [f"solid {name}"]
    for facet, normal in zip(tri, normals):
        lines.append("  facet normal {:.9e} {:.9e} {:.9e}".format(*normal))
        lines.append("    outer loop")
        for corner in facet:
            lines.append("      vertex {:.9e} {:.9e} {:.9e}".format(*corner))
        lines.append("    endloop")
        lines.append("  endfacet")
    lines.append(f"endsolid {name}")
    path.write_text("\n".join(lines) + "\n")
    return path


@dataclass(frozen=True)
class MeshDiagnostics:
    """Topological verdict on a welded triangle mesh."""

    n_vertices: int
    n_faces: int
    n_degenerate_faces: int
    n_boundary_edges: int
    n_nonmanifold_edges: int
    consistently_oriented: bool
    signed_volume: float

    @property
    def is_watertight(self):
        return self.n_boundary_edges == 0 and self.n_nonmanifold_edges == 0

    @property
    def outward_oriented(self):
        return self.signed_volume > 0.0

    def failure_reason(self):
        """Human-readable refusal text, or ``None`` when the mesh is usable."""
        if self.n_degenerate_faces:
            return (
                f"degenerate mesh: {self.n_degenerate_faces} zero-area "
                "facet(s); remove them with a topology-preserving repair"
            )
        if self.n_boundary_edges:
            return (f"not watertight: {self.n_boundary_edges} boundary edge(s) "
                    "bound only one facet (the surface has holes)")
        if self.n_nonmanifold_edges:
            return (f"non-manifold: {self.n_nonmanifold_edges} edge(s) bound "
                    "more than two facets")
        if not self.consistently_oriented:
            return ("inconsistent facet winding: shared edges are not traversed "
                    "in opposite directions")
        if self.signed_volume == 0.0:
            return "degenerate mesh: zero enclosed volume"
        return None


def diagnose_mesh(mesh):
    """Watertightness, manifoldness, winding consistency, and enclosed volume."""
    faces = np.asarray(mesh.faces, dtype=int)
    keep = mesh.face_areas > 0.0
    live = faces[keep]
    directed = np.concatenate((live[:, [0, 1]], live[:, [1, 2]], live[:, [2, 0]]))
    undirected = np.sort(directed, axis=1)
    _, counts = np.unique(undirected, axis=0, return_counts=True)
    _, directed_counts = np.unique(directed, axis=0, return_counts=True)
    return MeshDiagnostics(
        n_vertices=int(len(mesh.vertices)),
        n_faces=int(len(faces)),
        n_degenerate_faces=int(np.count_nonzero(~keep)),
        n_boundary_edges=int(np.count_nonzero(counts == 1)),
        n_nonmanifold_edges=int(np.count_nonzero(counts > 2)),
        consistently_oriented=bool(np.all(directed_counts == 1)
                                   and np.all(counts <= 2)),
        signed_volume=float(mesh.signed_volume))


def drop_degenerate_faces(mesh, *, relative_area_tolerance=0.0):
    """Remove only zero/near-zero facets and prove topology is unchanged.

    Some CAD exporters append point facets ``[v, v, v]``.  Ignoring them
    silently is unsafe because a zero-area facet can also be evidence of a
    collapsed or corrupted surface.  This helper removes the requested area
    class, compacts unused vertices, and accepts the repair only when the
    remaining surface is closed, manifold, consistently oriented, nonzero in
    volume, and volume-preserving to floating-point precision.

    Returns ``(clean_mesh, receipt)`` where ``receipt`` is a plain mapping that
    can be serialized beside the supplied STL.
    """
    tolerance = float(relative_area_tolerance)
    if not np.isfinite(tolerance) or tolerance < 0.0:
        raise ValueError("relative_area_tolerance must be finite and nonnegative")
    vertices = np.asarray(mesh.vertices, dtype=float)
    extent = float(np.max(np.ptp(vertices, axis=0)))
    area_limit = tolerance * extent ** 2
    areas = np.asarray(mesh.face_areas, dtype=float)
    removable = areas <= area_limit
    indices = np.flatnonzero(removable)
    if len(indices) == 0:
        return mesh, {
            "operation": "drop_degenerate_faces",
            "removed_face_indices": [],
            "removed_face_count": 0,
            "relative_area_tolerance": tolerance,
            "absolute_area_tolerance_file_units_squared": area_limit,
            "signed_volume_before": float(mesh.signed_volume),
            "signed_volume_after": float(mesh.signed_volume),
            "topology_preserved": True,
        }
    surviving_faces = np.asarray(mesh.faces, dtype=int)[~removable]
    if len(surviving_faces) == 0:
        raise ValueError("degenerate-face repair would remove the whole mesh")
    used, inverse = np.unique(surviving_faces.reshape(-1), return_inverse=True)
    compact_faces = inverse.reshape(-1, 3)
    normals = None
    if mesh.file_normals is not None:
        normals = np.asarray(mesh.file_normals, dtype=float)[~removable]
    clean = StlMesh(vertices[used], compact_faces, normals)
    diagnostics = diagnose_mesh(clean)
    reason = diagnostics.failure_reason()
    if reason is not None:
        raise ValueError(
            "degenerate-face removal did not preserve a usable surface -- "
            f"{reason}"
        )
    before = float(mesh.signed_volume)
    after = float(clean.signed_volume)
    scale = max(abs(before), abs(after), np.finfo(float).tiny)
    relative_volume_change = abs(after - before) / scale
    if relative_volume_change > 128.0 * np.finfo(float).eps:
        raise ValueError(
            "degenerate-face removal changed enclosed volume by "
            f"{relative_volume_change:.3e}"
        )
    return clean, {
        "operation": "drop_degenerate_faces",
        "removed_face_indices": indices.astype(int).tolist(),
        "removed_face_count": int(len(indices)),
        "relative_area_tolerance": tolerance,
        "absolute_area_tolerance_file_units_squared": area_limit,
        "signed_volume_before": before,
        "signed_volume_after": after,
        "relative_volume_change": relative_volume_change,
        "topology_preserved": True,
        "clean_diagnostics": {
            "n_vertices": diagnostics.n_vertices,
            "n_faces": diagnostics.n_faces,
            "n_degenerate_faces": diagnostics.n_degenerate_faces,
            "n_boundary_edges": diagnostics.n_boundary_edges,
            "n_nonmanifold_edges": diagnostics.n_nonmanifold_edges,
            "consistently_oriented": diagnostics.consistently_oriented,
            "signed_volume": diagnostics.signed_volume,
        },
    }


def _point_triangle_distance_and_solid_angle(points, tri):
    """Exact distances and summed solid angle from ``points`` to a triangle soup.

    Returns ``(distance, winding)`` where ``distance`` is the unsigned distance
    to the nearest triangle (Ericson's region-based closest point, vectorized)
    and ``winding`` is the Van Oosterom-Strackee solid-angle sum divided by
    4 pi -- exactly 1 inside a closed outward-oriented surface, 0 outside.
    """
    points = np.asarray(points, dtype=float)
    a, b, c = (np.asarray(tri[:, index], dtype=float) for index in range(3))
    ab = b - a
    ac = c - a
    bc = c - b
    pa = points[:, None, :] - a[None]
    d1 = np.einsum("fj,pfj->pf", ab, pa)
    d2 = np.einsum("fj,pfj->pf", ac, pa)
    pb = points[:, None, :] - b[None]
    d3 = np.einsum("fj,pfj->pf", ab, pb)
    d4 = np.einsum("fj,pfj->pf", ac, pb)
    pc = points[:, None, :] - c[None]
    d5 = np.einsum("fj,pfj->pf", ab, pc)
    d6 = np.einsum("fj,pfj->pf", ac, pc)
    va = d3 * d6 - d5 * d4
    vb = d5 * d2 - d1 * d6
    vc = d1 * d4 - d3 * d2
    with np.errstate(divide="ignore", invalid="ignore"):
        # Interior of the face (barycentric clamp).
        denominator = va + vb + vc
        safe = np.where(np.abs(denominator) > 0.0, denominator, 1.0)
        v_face = np.clip(vb / safe, 0.0, 1.0)
        w_face = np.clip(vc / safe, 0.0, 1.0)
        closest = (a[None] + v_face[..., None] * ab[None]
                   + w_face[..., None] * ac[None])
        # Edge AB.
        t_ab = np.clip(np.where(d1 - d3 != 0.0, d1 / np.where(d1 - d3 != 0.0, d1 - d3, 1.0),
                                0.0), 0.0, 1.0)
        edge_ab = a[None] + t_ab[..., None] * ab[None]
        # Edge AC.
        t_ac = np.clip(np.where(d2 - d6 != 0.0, d2 / np.where(d2 - d6 != 0.0, d2 - d6, 1.0),
                                0.0), 0.0, 1.0)
        edge_ac = a[None] + t_ac[..., None] * ac[None]
        # Edge BC.
        span = (d4 - d3) + (d5 - d6)
        t_bc = np.clip(np.where(span != 0.0, (d4 - d3) / np.where(span != 0.0, span, 1.0),
                                0.0), 0.0, 1.0)
        edge_bc = b[None] + t_bc[..., None] * bc[None]
    # Region select.  Written as an overwrite cascade in reverse priority, so
    # the surviving branch is the one Ericson's early-return order reaches
    # first: vertices beat edges, and edge AB beats AC beats BC.
    candidate = closest
    candidate = np.where(
        ((va <= 0.0) & ((d4 - d3) >= 0.0) & ((d5 - d6) >= 0.0))[..., None],
        edge_bc, candidate)
    candidate = np.where(
        ((vb <= 0.0) & (d2 >= 0.0) & (d6 <= 0.0))[..., None], edge_ac, candidate)
    candidate = np.where(
        ((vc <= 0.0) & (d1 >= 0.0) & (d3 <= 0.0))[..., None], edge_ab, candidate)
    candidate = np.where(
        ((d6 >= 0.0) & (d5 <= d6))[..., None], c[None], candidate)
    candidate = np.where(
        ((d3 >= 0.0) & (d4 <= d3))[..., None], b[None], candidate)
    candidate = np.where(
        ((d1 <= 0.0) & (d2 <= 0.0))[..., None], a[None], candidate)
    distance = np.sqrt(np.sum((points[:, None, :] - candidate) ** 2, axis=2))
    # Solid angle of each triangle as seen from each point.
    la = np.linalg.norm(pa, axis=2)
    lb = np.linalg.norm(pb, axis=2)
    lc = np.linalg.norm(pc, axis=2)
    # Van Oosterom-Strackee is written in the vectors (A - P); pa..pc above are
    # (P - A), whose triple product is the negation (the pairwise dot products
    # and norms in the denominator are unchanged by the flip).
    numerator = -np.einsum("pfj,pfj->pf", pa, np.cross(pb, pc))
    denominator = (la * lb * lc
                   + np.einsum("pfj,pfj->pf", pa, pb) * lc
                   + np.einsum("pfj,pfj->pf", pb, pc) * la
                   + np.einsum("pfj,pfj->pf", pc, pa) * lb)
    omega = 2.0 * np.arctan2(numerator, denominator)
    return (np.min(distance, axis=1), np.sum(omega, axis=1) / (4.0 * np.pi))


@dataclass(frozen=True)
class SdfRasterReport:
    """Receipts for one STL -> signed-distance rasterization."""

    shape: tuple
    dx: float
    origin: tuple
    diagnostics: MeshDiagnostics
    solid_region: str
    solid_fraction: float
    min_feature_cells: float
    max_winding_residual: float

    def resolution_warning(self, required_cells=6.0):
        """Text warning when the narrowest resolved band is under-resolved."""
        if self.min_feature_cells >= required_cells:
            return None
        return (f"narrowest resolved band is {self.min_feature_cells:.1f} cells "
                f"(< {required_cells:.0f}); refine dx or the feature width will "
                "not converge")


def _grid_points(shape, dx, origin):
    axes = [origin[axis] + np.arange(shape[axis]) * dx for axis in range(3)]
    grid = np.meshgrid(*axes, indexing="ij")
    return np.stack([item.reshape(-1) for item in grid], axis=1)


def rasterize_signed_distance(
        mesh, *, dx, shape=None, origin=None, padding_cells=2.0,
        solid_region="interior", solid_ceiling=None, solid_floor=None,
        chunk_pairs=4_000_000, require_watertight=True):
    """Rasterize a watertight STL to a petch signed distance field.

    Returns ``(phi, report)`` with ``phi > 0`` inside the declared solid, in the
    same length units as the STL and ``dx``.  ``solid_region="interior"`` treats
    the STL body as solid; ``"exterior"`` treats it as the void (the hole), so
    the solid fills the remainder of the grid box.  ``solid_ceiling`` /
    ``solid_floor`` intersect the result with the half-spaces ``z <= ceiling`` /
    ``z >= floor``, which is how an open feature mouth is declared.
    """
    if solid_region not in {"interior", "exterior"}:
        raise ValueError("solid_region must be 'interior' or 'exterior'")
    diagnostics = diagnose_mesh(mesh)
    reason = diagnostics.failure_reason()
    if require_watertight and reason is not None:
        raise ValueError(f"STL cannot be rasterized -- {reason}")
    oriented = mesh if diagnostics.signed_volume > 0.0 else mesh.flipped()
    lower, upper = (np.asarray(bound, dtype=float) for bound in oriented.bounds)
    dx = float(dx)
    if not np.isfinite(dx) or dx <= 0.0:
        raise ValueError("dx must be finite and positive")
    if origin is None:
        # Stagger the default grid by half a cell.  CAD bodies are overwhelmingly
        # axis-aligned on round coordinates, so an unstaggered origin puts whole
        # facet planes exactly on grid nodes: phi is then exactly zero there,
        # which is neither solid nor gas, and marching cubes returns a surface
        # with unmatched edges.  Half-cell stagger removes the coincidence for
        # any facet plane separated from the bounds by a multiple of dx.
        origin = tuple(lower - (padding_cells + 0.5) * dx)
    origin = tuple(float(value) for value in origin)
    if shape is None:
        span = (upper + padding_cells * dx) - np.asarray(origin)
        shape = tuple(int(np.ceil(value / dx)) + 1 for value in span)
    shape = tuple(int(value) for value in shape)
    if len(shape) != 3 or min(shape) < 2:
        raise ValueError("grid shape must be 3-D with at least two cells per axis")

    points = _grid_points(shape, dx, origin)
    tri = np.asarray(oriented.triangles, dtype=float)
    keep = oriented.face_areas > 0.0
    tri = tri[keep]
    if len(tri) == 0:
        raise ValueError("STL has no non-degenerate facets")
    chunk = max(1, int(chunk_pairs // len(tri)))
    distance = np.empty(len(points))
    winding = np.empty(len(points))
    for start in range(0, len(points), chunk):
        stop = min(start + chunk, len(points))
        block_distance, block_winding = _point_triangle_distance_and_solid_angle(
            points[start:stop], tri)
        distance[start:stop] = block_distance
        winding[start:stop] = block_winding
    inside = winding > 0.5
    residual = float(np.max(np.minimum(np.abs(winding), np.abs(winding - 1.0))))
    signed = np.where(inside, distance, -distance)
    if solid_region == "exterior":
        signed = -signed
    phi = signed.reshape(shape)
    if solid_ceiling is not None:
        z = origin[2] + np.arange(shape[2]) * dx
        phi = np.minimum(phi, float(solid_ceiling) - z[None, None, :])
    if solid_floor is not None:
        z = origin[2] + np.arange(shape[2]) * dx
        phi = np.minimum(phi, z[None, None, :] - float(solid_floor))
    solid_fraction = float(np.count_nonzero(phi > 0.0) / phi.size)
    report = SdfRasterReport(
        shape=shape, dx=dx, origin=origin, diagnostics=diagnostics,
        solid_region=solid_region, solid_fraction=solid_fraction,
        min_feature_cells=_min_feature_cells(phi, dx),
        max_winding_residual=residual)
    return phi, report


def _min_feature_cells(phi, dx):
    """Narrowest resolved solid or gas run along any axis, in cells.

    A run of one cell means the band exists only as a single sample and is not
    representable by the level set; this is the receipt that lets a caller
    refuse an under-resolved import instead of etching a different structure.
    """
    best = np.inf
    solid = phi > 0.0
    for axis in range(3):
        lines = np.moveaxis(solid, axis, 0).reshape(solid.shape[axis], -1).T
        length = lines.shape[1]
        for values in (lines, ~lines):
            padded = np.zeros((values.shape[0], length + 2), dtype=np.int8)
            padded[:, 1:-1] = values
            changes = np.diff(padded, axis=1)
            starts = np.argwhere(changes == 1)
            stops = np.argwhere(changes == -1)
            if not len(starts):
                continue
            # argwhere is row-major, so run i on a line pairs start i with stop i.
            # A run touching either end of the line is truncated by the domain,
            # not by geometry -- counting it would report the padding width as a
            # feature size.
            interior = (starts[:, 1] > 0) & (stops[:, 1] < length)
            if np.any(interior):
                best = min(best, float(np.min(
                    (stops[interior, 1] - starts[interior, 1]))))
    return float(best) if np.isfinite(best) else float(min(phi.shape))


@dataclass(frozen=True)
class AxisymmetryReport:
    """Measured generator profile ``r(z)`` plus its axisymmetry deviation."""

    z: np.ndarray
    r: np.ndarray
    level_deviation: np.ndarray
    axis_point: tuple
    axis_direction: tuple
    max_deviation: float
    relative_deviation: float
    facet_bound: float
    is_axisymmetric: bool

    def profile(self, **kwargs):
        return to_axisymmetric_profile(self, **kwargs)


def _axis_frame(vertices, axis_point, axis_direction):
    direction = np.asarray(axis_direction, dtype=float)
    norm = float(np.linalg.norm(direction))
    if not np.isfinite(norm) or norm == 0.0:
        raise ValueError("axis_direction must be a nonzero vector")
    direction = direction / norm
    if axis_point is None:
        lower = vertices.min(axis=0)
        upper = vertices.max(axis=0)
        axis_point = 0.5 * (lower + upper)
    axis_point = np.asarray(axis_point, dtype=float)
    offset = vertices - axis_point
    radial = offset - (offset @ direction)[:, None] * direction[None, :]
    # Height is the absolute projection onto the axis direction (for the default
    # z axis this is simply the z coordinate), so the extracted profile shares
    # the STL's own coordinates instead of being shifted to the axis point.
    height = vertices @ direction
    return axis_point, direction, height, np.linalg.norm(radial, axis=1)


def _facet_bound(mesh, axis_point, direction):
    """Radius error implied by the azimuthal tessellation of this mesh.

    Estimated from the largest azimuthal step between adjacent welded vertices
    at a common height, i.e. ``R (1 - cos(dphi / 2))`` for the coarsest band.
    """
    vertices = np.asarray(mesh.vertices, dtype=float)
    offset = vertices - axis_point
    height = offset @ direction
    radial = offset - height[:, None] * direction[None, :]
    radius = np.linalg.norm(radial, axis=1)
    live = radius > 0.0
    if np.count_nonzero(live) < 3:
        return 0.0
    basis = np.eye(3)[np.argmin(np.abs(direction))]
    e1 = np.cross(direction, basis)
    e1 = e1 / np.linalg.norm(e1)
    e2 = np.cross(direction, e1)
    angle = np.arctan2(radial[live] @ e2, radial[live] @ e1)
    unique = np.unique(np.round(angle, 9))
    if unique.size < 2:
        return 0.0
    steps = np.diff(np.sort(unique))
    largest = float(max(steps.max(), 2.0 * np.pi - float(np.sum(steps))))
    return float(np.max(radius) * (1.0 - np.cos(0.5 * largest)))


def extract_axisymmetric_profile(
        mesh, *, axis_point=None, axis_direction=(0.0, 0.0, 1.0), n_levels=48,
        relative_tolerance=None, tolerance_safety=4.0,
        max_relative_deviation=1e-2):
    """Slice a body of revolution into a generator profile ``r(z)``.

    Each level is cut by intersecting every crossing mesh edge with the plane,
    so the extracted radius is exact for a piecewise-linear body of revolution
    and the reported ``max_deviation`` is the true out-of-roundness of the file
    (tessellation faceting included, not assumed away).

    Acceptance requires the measured out-of-roundness to be consistent with the
    mesh's own faceting bound *and* below ``max_relative_deviation``.  Both are
    needed: the faceting bound alone accepts a coarse polygon (a square prism
    has only four azimuthal samples, so its bound is ~29% of the radius and
    would excuse anything), while the cap alone would reject a legitimately
    coarse tessellation without saying why.  Pinning ``relative_tolerance``
    overrides both with a declared number.
    """
    vertices = np.asarray(mesh.vertices, dtype=float)
    axis_point, direction, height, _ = _axis_frame(
        vertices, axis_point, axis_direction)
    faces = np.asarray(mesh.faces, dtype=int)
    edges = np.unique(np.sort(np.concatenate(
        (faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]])), axis=1), axis=0)
    low = float(height.min())
    high = float(height.max())
    if not np.isfinite(low) or high <= low:
        raise ValueError("degenerate extent along the declared axis")
    n_levels = int(n_levels)
    if n_levels < 2:
        raise ValueError("n_levels must be at least 2")
    inset = 0.5 * (high - low) / n_levels
    levels = np.linspace(low + inset, high - inset, n_levels)

    z_a = height[edges[:, 0]]
    z_b = height[edges[:, 1]]
    v_a = vertices[edges[:, 0]]
    v_b = vertices[edges[:, 1]]
    radii = []
    deviation = []
    kept_levels = []
    for level in levels:
        crossing = ((z_a - level) * (z_b - level) < 0.0)
        if not np.any(crossing):
            continue
        span = z_b[crossing] - z_a[crossing]
        fraction = (level - z_a[crossing]) / span
        point = v_a[crossing] + fraction[:, None] * (v_b[crossing] - v_a[crossing])
        offset = point - axis_point
        along = offset @ direction
        radial = np.linalg.norm(offset - along[:, None] * direction[None, :], axis=1)
        mean_radius = float(np.mean(radial))
        radii.append(mean_radius)
        deviation.append(float(np.max(np.abs(radial - mean_radius))))
        kept_levels.append(float(level))
    if len(kept_levels) < 2:
        raise ValueError("fewer than two levels intersect the mesh")
    radii = np.asarray(radii)
    deviation = np.asarray(deviation)
    scale = float(np.mean(radii))
    max_deviation = float(np.max(deviation))
    relative = max_deviation / scale if scale > 0.0 else np.inf
    bound = _facet_bound(mesh, axis_point, direction)
    if relative_tolerance is not None:
        threshold = float(relative_tolerance)
    else:
        faceting = tolerance_safety * (bound / scale if scale > 0.0 else 0.0)
        threshold = min(float(max_relative_deviation), faceting)
    return AxisymmetryReport(
        z=_readonly(kept_levels), r=_readonly(radii),
        level_deviation=_readonly(deviation),
        axis_point=tuple(float(value) for value in axis_point),
        axis_direction=tuple(float(value) for value in direction),
        max_deviation=max_deviation, relative_deviation=float(relative),
        facet_bound=float(bound),
        is_axisymmetric=bool(relative <= max(threshold, 1e-12)))


def to_axisymmetric_profile(report, *, require_axisymmetric=True):
    """Convert a report into the validated operator's generator profile.

    The last vertex is the opening rim, matching
    :class:`petch.axisymmetric_exchange_3d.AxisymmetricProfile`.
    """
    if require_axisymmetric and not report.is_axisymmetric:
        raise ValueError(
            "geometry is not axisymmetric: relative deviation "
            f"{report.relative_deviation:.3e} exceeds the mesh faceting bound "
            f"{report.facet_bound:.3e} -- use the 3-D level-set path")
    return AxisymmetricProfile(np.asarray(report.z), np.asarray(report.r))


def revolved_stl_mesh(z, r, *, n_theta=64, axis_point=(0.0, 0.0), close_bottom=True,
                      close_top=True):
    """Build a closed body-of-revolution mesh from a generator profile.

    The inverse of :func:`extract_axisymmetric_profile`: it turns an
    ``(z, r)`` generator (a hole or pillar profile) into a watertight,
    outward-wound STL mesh, so a validated axisymmetric profile can be exported
    to the 3-D path or written to file.
    """
    z = np.asarray(z, dtype=float)
    r = np.asarray(r, dtype=float)
    if z.shape != r.shape or z.ndim != 1 or len(z) < 2:
        raise ValueError("z and r must be matching 1-D generator arrays")
    if np.any(np.diff(z) <= 0.0) or np.any(r < 0.0):
        raise ValueError("z must strictly increase and r must be nonnegative")
    n_theta = int(n_theta)
    if n_theta < 3:
        raise ValueError("n_theta must be at least 3")
    theta = np.linspace(0.0, 2.0 * np.pi, n_theta, endpoint=False)
    cos = np.cos(theta)
    sin = np.sin(theta)
    ring = np.stack([
        np.stack([axis_point[0] + radius * cos,
                  axis_point[1] + radius * sin,
                  np.full(n_theta, level)], axis=1)
        for level, radius in zip(z, r)])
    triangles = []
    for index in range(len(z) - 1):
        lower = ring[index]
        upper = ring[index + 1]
        for j in range(n_theta):
            k = (j + 1) % n_theta
            triangles.append([lower[j], lower[k], upper[k]])
            triangles.append([lower[j], upper[k], upper[j]])
    if close_bottom and r[0] > 0.0:
        center = np.array([axis_point[0], axis_point[1], z[0]])
        for j in range(n_theta):
            k = (j + 1) % n_theta
            triangles.append([center, ring[0][k], ring[0][j]])
    if close_top and r[-1] > 0.0:
        center = np.array([axis_point[0], axis_point[1], z[-1]])
        for j in range(n_theta):
            k = (j + 1) % n_theta
            triangles.append([center, ring[-1][j], ring[-1][k]])
    corners = np.asarray(triangles, dtype=float)
    areas = 0.5 * np.linalg.norm(np.cross(
        corners[:, 1] - corners[:, 0], corners[:, 2] - corners[:, 0]), axis=1)
    corners = corners[areas > 0.0]
    vertices, faces = _weld(corners)
    return StlMesh(vertices, faces)


def assign_materials_by_z(phi, *, dx, origin=(0.0, 0.0, 0.0), layers):
    """Split a solid level set into stacked materials by height.

    ``layers`` is a sequence of ``(z_low, z_high, material_id)`` in the same
    units as ``dx``, contiguous and covering the grid extent; ids must be
    positive (zero is gas).  Each material level set is the CSG intersection of
    the solid with its slab, so their pointwise maximum reproduces ``phi``'s
    sign exactly -- the invariant :class:`FeatureGeometry3D` enforces.
    """
    phi = np.asarray(phi, dtype=float)
    if phi.ndim != 3:
        raise ValueError("phi must be a 3-D field")
    entries = [(float(low), float(high), int(identifier))
               for low, high, identifier in layers]
    if not entries:
        raise ValueError("at least one material layer is required")
    if any(identifier <= 0 for _, _, identifier in entries):
        raise ValueError("material ids must be positive (0 is gas)")
    if len({identifier for _, _, identifier in entries}) != len(entries):
        raise ValueError("material ids must be unique")
    entries.sort()
    z = origin[2] + np.arange(phi.shape[2]) * dx
    for (_, high, _), (low, _, _) in zip(entries, entries[1:]):
        if abs(high - low) > 1e-12 * max(1.0, abs(high)):
            raise ValueError("material layers must be contiguous in z")
    if entries[0][0] > z[0] + 1e-12 or entries[-1][1] < z[-1] - 1e-12:
        raise ValueError("material layers must cover the grid extent in z")
    levelsets = {}
    for low, high, identifier in entries:
        slab = np.minimum(z - low, high - z)[None, None, :]
        levelsets[identifier] = np.minimum(phi, np.broadcast_to(slab, phi.shape))
    stacked = np.stack([levelsets[identifier] for _, _, identifier in entries])
    ids = np.asarray([identifier for _, _, identifier in entries])
    material = np.where(phi > 0.0, ids[np.argmax(stacked, axis=0)], 0)
    return material.astype(int), levelsets


def build_feature_geometry_from_stl(
        source, *, dx, mesh_length_unit_m, layers=None, material_id=1,
        shape=None, origin=None, padding_cells=2.0, solid_region="interior",
        solid_ceiling=None, solid_floor=None, require_watertight=True,
        chunk_pairs=4_000_000):
    """STL (path, bytes, or :class:`StlMesh`) -> engine-ready feature geometry.

    Returns ``(geometry, report)``.  With ``layers`` the solid is split into
    stacked materials (mask over substrate); without it the whole solid carries
    ``material_id``.  Lengths stay in STL units; ``mesh_length_unit_m`` declares
    what one unit is in metres, exactly as elsewhere in the engine.
    """
    mesh = source if isinstance(source, StlMesh) else read_stl(source)
    phi, report = rasterize_signed_distance(
        mesh, dx=dx, shape=shape, origin=origin, padding_cells=padding_cells,
        solid_region=solid_region, solid_ceiling=solid_ceiling,
        solid_floor=solid_floor, require_watertight=require_watertight,
        chunk_pairs=chunk_pairs)
    if layers is None:
        material = np.where(phi > 0.0, int(material_id), 0)
        levelsets = {int(material_id): phi}
    else:
        material, levelsets = assign_materials_by_z(
            phi, dx=report.dx, origin=report.origin, layers=layers)
    # Engine mesh coordinates are index * dx (FeatureGeometry3D.coordinate_arrays
    # ignores the origin); the STL-frame origin is carried as the SI offset of
    # index (0, 0, 0), which is what the transport reference-plane check reads.
    geometry = FeatureGeometry3D(
        phi, material, report.dx, float(mesh_length_unit_m),
        tuple(value * float(mesh_length_unit_m) for value in report.origin),
        material_levelsets=levelsets)
    return geometry, report
