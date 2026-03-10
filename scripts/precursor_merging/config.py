"""
Matching Configuration

Defines configuration classes and enums for precursor matching.
Moved from precursor_matching.py for modular organization.
"""

from dataclasses import dataclass
from enum import Enum
from typing import Optional


class MatchTier(Enum):
    """Match quality tiers.

    Higher tiers indicate higher confidence matches.
    Values are preserved for backwards compatibility.
    """
    SEQUENCE_EXACT = 1      # Exact sequence + charge match
    SEQUENCE_IL_NORM = 2    # I/L normalized sequence + charge match
    COORDINATE_FULL = 3     # m/z + charge + RT + IM match
    COORDINATE_PARTIAL = 4  # m/z + charge match (missing RT/IM)
    NO_MATCH = 5            # No match found


@dataclass
class MatchConfig:
    """Configuration for precursor matching.

    Default values are calibrated for timsTOF data:
    - 20 ppm m/z tolerance: Standard for high-resolution MS
    - 5s RT tolerance: Tight window within same run
    - 0.05 1/K0 IM tolerance: Typical for timsTOF precision
    """
    # m/z tolerance in ppm
    mz_tol_ppm: float = 20.0

    # RT tolerance in seconds
    # Tightened to 0.01s (2026-03-10): Vadim Demichev confirmed DIA-NN DDA RT
    # is frame-exact for timsTOF ("can even go to 0.01s"). 10ms is sub-frame
    # precision (~100ms frame rate), effectively an exact match.
    # History: 5.0s (caused 78%+ stranded) → 0.5s (worked) → 0.01s (frame-exact).
    rt_tol_sec: float = 0.01

    # Ion mobility tolerance (1/K0)
    im_tol: float = 0.05

    # Charge must match exactly
    charge_must_match: bool = True

    # Whether to normalize I/L in sequences for matching
    normalize_il: bool = True

    def __str__(self) -> str:
        return (
            f"MatchConfig(mz_tol={self.mz_tol_ppm}ppm, "
            f"rt_tol={self.rt_tol_sec}s, im_tol={self.im_tol})"
        )


@dataclass
class MatchResult:
    """Result of a precursor match attempt."""
    matched: bool
    tier: MatchTier
    score: float  # Lower is better

    # Differences (for diagnostics)
    mz_diff_ppm: Optional[float] = None
    rt_diff_sec: Optional[float] = None
    im_diff: Optional[float] = None
    sequence_match: bool = False
    charge_match: bool = False

    @classmethod
    def no_match(cls) -> 'MatchResult':
        """Factory for no-match result."""
        return cls(
            matched=False,
            tier=MatchTier.NO_MATCH,
            score=float('inf'),
        )

    @classmethod
    def sequence_match(
        cls,
        normalize_il: bool = True,
        mz_diff_ppm: Optional[float] = None,
        rt_diff_sec: Optional[float] = None,
        im_diff: Optional[float] = None,
    ) -> 'MatchResult':
        """Factory for sequence match result."""
        return cls(
            matched=True,
            tier=MatchTier.SEQUENCE_IL_NORM if normalize_il else MatchTier.SEQUENCE_EXACT,
            score=0.0,
            mz_diff_ppm=mz_diff_ppm,
            rt_diff_sec=rt_diff_sec,
            im_diff=im_diff,
            sequence_match=True,
            charge_match=True,
        )

    @classmethod
    def coordinate_match(
        cls,
        score: float,
        full: bool,
        mz_diff_ppm: float,
        rt_diff_sec: Optional[float] = None,
        im_diff: Optional[float] = None,
    ) -> 'MatchResult':
        """Factory for coordinate match result."""
        return cls(
            matched=True,
            tier=MatchTier.COORDINATE_FULL if full else MatchTier.COORDINATE_PARTIAL,
            score=score,
            mz_diff_ppm=mz_diff_ppm,
            rt_diff_sec=rt_diff_sec,
            im_diff=im_diff,
            sequence_match=False,
            charge_match=True,
        )
