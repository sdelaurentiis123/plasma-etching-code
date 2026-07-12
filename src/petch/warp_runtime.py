"""Shared Warp runtime setup that changes cache storage only, never device or numerics."""
from __future__ import annotations

import os
from pathlib import Path
import tempfile


_CONFIGURED = set()


def ensure_writable_warp_cache(wp):
    """Select a process-local cache when Warp's configured user cache is not writable."""
    key = id(wp)
    if key in _CONFIGURED:
        return Path(wp.config.kernel_cache_dir)
    cache = Path(wp.config.kernel_cache_dir)
    try:
        cache.mkdir(parents=True, exist_ok=True)
        probe = Path(tempfile.mkdtemp(prefix="petch-probe-", dir=cache))
        probe.rmdir()
    except OSError:
        fallback = Path(tempfile.gettempdir()) / "petch-warp-cache" / wp.config.version
        fallback.mkdir(parents=True, exist_ok=True)
        wp.config.kernel_cache_dir = os.fspath(fallback)
        cache = fallback
    _CONFIGURED.add(key)
    return cache
