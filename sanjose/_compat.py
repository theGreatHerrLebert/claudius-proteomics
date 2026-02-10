"""Optional dependency checks."""

HAS_ZSTD = False
try:
    import zstandard as zstd
    HAS_ZSTD = True
except ImportError:
    zstd = None


def get_zstd_decompressor():
    """Return a zstd decompressor, raising ImportError if unavailable."""
    if not HAS_ZSTD:
        raise ImportError(
            "zstandard is required for reading zstd-compressed blobs. "
            "Install with: pip install zstandard"
        )
    return zstd.ZstdDecompressor()
