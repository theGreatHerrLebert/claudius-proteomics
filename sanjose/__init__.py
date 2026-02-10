"""San José: Python query API for timsTOF precursor data.

Usage:
    from sanjose import Dataset

    ds = Dataset("data/merged/PXD019086")
    df = ds.query(sequence="YRPGTVALR", charge=2)

    for p in ds.iter(min_engines=3, limit=10):
        sig = p.signal
        print(p.raw_file, sig.fragment_spectrum.n_peaks)
"""

from .dataset import Dataset
from .store import PrecursorStore
from .blob import BlobReader
from .models import Precursor, FragmentSpectrum, MS1Signal, RawPointCloud, RawSignal

__all__ = [
    'Dataset',
    'PrecursorStore',
    'BlobReader',
    'Precursor',
    'FragmentSpectrum',
    'MS1Signal',
    'RawPointCloud',
    'RawSignal',
]
