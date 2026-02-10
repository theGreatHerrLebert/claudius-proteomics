"""Data models for San José precursor data."""

from dataclasses import dataclass, field
from typing import Optional, Dict, Any, TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from .blob import BlobReader


@dataclass
class FragmentSpectrum:
    """Merged MS2 fragment spectrum."""
    mz: np.ndarray
    intensity: np.ndarray
    mobility: Optional[np.ndarray] = None
    scan: Optional[np.ndarray] = None

    @property
    def n_peaks(self) -> int:
        return len(self.mz)

    def __repr__(self) -> str:
        return f"<FragmentSpectrum: {self.n_peaks} peaks>"


@dataclass
class MS1Signal:
    """1D projections of MS1 precursor signal."""
    xic_rt: np.ndarray
    xic_intensity: np.ndarray
    mobilogram_im: np.ndarray
    mobilogram_intensity: np.ndarray
    isotope_intensity: np.ndarray

    def __repr__(self) -> str:
        return (
            f"<MS1Signal: XIC={len(self.xic_rt)}pts, "
            f"mobilogram={len(self.mobilogram_im)}pts, "
            f"isotopes={len(self.isotope_intensity)}>"
        )


@dataclass
class RawPointCloud:
    """Raw 4D MS1 point cloud (all peaks in RT/IM window)."""
    rt: np.ndarray
    mz: np.ndarray
    mobility: np.ndarray
    intensity: np.ndarray

    @property
    def n_points(self) -> int:
        return len(self.rt)

    def __repr__(self) -> str:
        return f"<RawPointCloud: {self.n_points} points>"


@dataclass
class RawSignal:
    """Complete raw signal for a precursor (loaded from blob)."""
    fragment_spectrum: FragmentSpectrum
    ms1_signal: MS1Signal
    raw_point_cloud: RawPointCloud
    metadata: Dict[str, Any] = field(default_factory=dict)

    def __repr__(self) -> str:
        return (
            f"<RawSignal: {self.fragment_spectrum.n_peaks} fragments, "
            f"{self.raw_point_cloud.n_points} raw points>"
        )


@dataclass
class Precursor:
    """A single precursor observation with lazy signal loading.

    Scalar fields are always available. The raw 4D signal is loaded
    on first access of the `.signal` property from the blob file.
    """
    # Identity
    precursor_id: int
    raw_file: str

    # Measured properties
    mono_mz: float
    charge: int
    rt_seconds: float
    mobility: float

    # Engine identifications
    n_engines: int = 0
    consensus_peptide: Optional[str] = None
    fragpipe_peptide: Optional[str] = None
    fragpipe_modified: Optional[str] = None
    diann_peptide: Optional[str] = None
    diann_modified: Optional[str] = None
    sage_peptide: Optional[str] = None
    sage_modified: Optional[str] = None

    # Quality metrics
    is_high_quality: Optional[bool] = None
    ms1_rt_r2: Optional[float] = None
    ms1_im_r2: Optional[float] = None
    isotope_cosim: Optional[float] = None

    # Engine scores
    fragpipe_probability: Optional[float] = None
    diann_qvalue: Optional[float] = None
    sage_qvalue: Optional[float] = None

    # Blob reference (internal, for lazy loading)
    _blob_offset: Optional[int] = field(default=None, repr=False)
    _blob_size: Optional[int] = field(default=None, repr=False)
    _blob_reader: Optional['BlobReader'] = field(default=None, repr=False)
    _signal_cache: Optional[RawSignal] = field(default=None, repr=False)

    @property
    def signal(self) -> Optional[RawSignal]:
        """Load raw 4D signal from blob on first access."""
        if self._signal_cache is not None:
            return self._signal_cache

        if self._blob_reader is None or self._blob_offset is None or self._blob_size is None:
            return None

        self._signal_cache = self._blob_reader.read(
            self.raw_file, self._blob_offset, self._blob_size
        )
        return self._signal_cache

    @property
    def has_blob(self) -> bool:
        return self._blob_offset is not None and self._blob_size is not None

    def __repr__(self) -> str:
        seq = self.consensus_peptide or "unidentified"
        blob = ", blob" if self.has_blob else ""
        return (
            f"<Precursor {self.raw_file}:{self.precursor_id} "
            f"{seq}/{self.charge}+ "
            f"mz={self.mono_mz:.4f} "
            f"engines={self.n_engines}{blob}>"
        )
