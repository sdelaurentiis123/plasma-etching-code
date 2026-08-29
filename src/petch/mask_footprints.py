"""Analytic periodic mask-footprint fields for layered 3-D feature geometry."""
from __future__ import annotations

import numpy as np


def _centered_grid(*, cell_width, cell_length, dx):
    values = np.asarray([cell_width, cell_length, dx], dtype=float)
    if np.any(~np.isfinite(values)) or np.any(values <= 0.0):
        raise ValueError("mask-footprint cell dimensions must be positive")
    shape = tuple(
        max(3, int(round(length / dx)) + 1)
        for length in (cell_width, cell_length)
    )
    x = np.arange(shape[0]) * dx - 0.5 * float(cell_width)
    y = np.arange(shape[1]) * dx - 0.5 * float(cell_length)
    return np.meshgrid(x, y, indexing="ij")


def centered_rectangle_footprint_levelset(
        *, cell_width, cell_length, dx, rectangle_width,
        rectangle_length):
    """Positive-inside field for a centered rectangular mask island."""
    values = np.asarray([rectangle_width, rectangle_length], dtype=float)
    if (
        np.any(~np.isfinite(values))
        or np.any(values <= 0.0)
        or rectangle_width >= cell_width
        or rectangle_length >= cell_length
    ):
        raise ValueError("rectangle must fit strictly inside the periodic cell")
    X, Y = _centered_grid(
        cell_width=cell_width, cell_length=cell_length, dx=dx)
    return np.minimum(
        0.5 * float(rectangle_width) - np.abs(X),
        0.5 * float(rectangle_length) - np.abs(Y),
    )


def centered_square_footprint_levelset(*, pitch, dx, square_width):
    """Positive-inside field for a centered square mask island."""
    return centered_rectangle_footprint_levelset(
        cell_width=pitch,
        cell_length=pitch,
        dx=dx,
        rectangle_width=square_width,
        rectangle_length=square_width,
    )


def centered_cross_footprint_levelset(
        *, pitch, dx, outer_width, arm_width):
    """Positive-inside field for the union of orthogonal rectangles."""
    values = np.asarray([outer_width, arm_width], dtype=float)
    if (
        np.any(~np.isfinite(values))
        or np.any(values <= 0.0)
        or arm_width >= outer_width
        or outer_width >= pitch
    ):
        raise ValueError("cross dimensions must satisfy 0 < arm < outer < pitch")
    horizontal = centered_rectangle_footprint_levelset(
        cell_width=pitch,
        cell_length=pitch,
        dx=dx,
        rectangle_width=outer_width,
        rectangle_length=arm_width,
    )
    vertical = centered_rectangle_footprint_levelset(
        cell_width=pitch,
        cell_length=pitch,
        dx=dx,
        rectangle_width=arm_width,
        rectangle_length=outer_width,
    )
    return np.maximum(horizontal, vertical)


def centered_inverse_square_hole_footprint_levelset(
        *, pitch, dx, opening_width):
    """Positive-in-mask field for a blanket mask with a centered square hole."""
    if (
        not np.isfinite(opening_width)
        or opening_width <= 0.0
        or opening_width >= pitch
    ):
        raise ValueError("opening must fit strictly inside the periodic cell")
    X, Y = _centered_grid(cell_width=pitch, cell_length=pitch, dx=dx)
    half = 0.5 * float(opening_width)
    return np.maximum(np.abs(X) - half, np.abs(Y) - half)
