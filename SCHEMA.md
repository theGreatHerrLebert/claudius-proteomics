# Schema — Claudius timsTOF DDA-PASEF PSM corpus

Field-by-field reference for both configs. All files are zstd-compressed parquet.
The `split` column (`train` / `validation` / `test`) is present in every row
(peptide-hash on `sequence_normalized`, seed 0).

**Reading conventions**
- *null* never means zero — it means the value was not available (e.g. an engine
  did not identify this precursor, or a fit failed). Distinguish with the
  explicit status/flag columns.
- The inclusion floor is a **union of per-engine reported q ≤ 0.01** (rank-1
  target), **not** a corpus-level FDR. Use the exposed scores to filter.

---

## `tier1_psms` — one row per accepted precursor identification

Grain: `(raw_file, precursor_id)`.

### Identity
| field | type | notes |
|---|---|---|
| `accession` | str | PRIDE PXD accession |
| `raw_file` | str | Bruker `.d` run name |
| `precursor_id` | int64 | precursor id within the run |
| `sage_psm_id` | float64 | Sage PSM id (null if FragPipe-only) |
| `sequence` | str | plain amino-acid sequence |
| `modified_sequence` | str | canonical engine's modified sequence (ProForma-style `[UNIMOD:N]`) |
| `sequence_normalized` | str | unmodified, lower-case — the cross-engine + split key |
| `charge` | int64 | precursor charge |
| `protein` | str | protein accession(s) |

### Engine status & confidence
| field | type | notes |
|---|---|---|
| `n_engines` | int64 | 1 or 2 engines supporting the canonical peptidoform+charge |
| `peptidoform_conflict` | bool | both engines passed the floor but assigned different peptidoforms |
| `sage_assignment_passes_floor` / `fragpipe_assignment_passes_floor` | bool | per-engine floor pass |
| `sage_qvalue`, `sage_pep`, `sage_hyperscore`, `sage_cosine` | float64 | Sage scores (null if not identified by Sage) |
| `fragpipe_qvalue`, `fragpipe_pep`, `fragpipe_probability`, `fragpipe_hyperscore` | float64 | FragPipe scores (null if not identified) |

### Apex estimates
| field | type | notes |
|---|---|---|
| `mz` | float64 | precursor m/z |
| `rt_seconds` | float64 | raw RT apex (s, per-run scale) |
| `rt_aligned` | float64 | **Sage cross-run-aligned RT (~[0,1]) — prefer this for RT modeling** (null if FragPipe-only) |
| `mobility` | float64 | ion mobility, 1/K0 |
| `sage_rt` / `sage_mobility` / `fragpipe_rt` / `fragpipe_mobility` | float64 | per-engine apex |

### Collision energy (NOVEL — raw volts)
| field | type | notes |
|---|---|---|
| `collision_energy_mean_v` | float64 | mean CE over the precursor's PASEF events (V) |
| `collision_energies_v` | list&lt;float64&gt; | per-PASEF-event CE (V); divide by 100 for the common normalized convention |

### RT / IM peak-shape labels
| field | type | notes |
|---|---|---|
| `ms1_rt_apex/fwhm/sigma/skew/r2` | float64 | RT peak fit (s); `sigma` is a Gaussian-LSQ width (over-widens for tailing peaks) |
| `ms1_rt_snr` | float64 | RT trace SNR (`trace_snr_v1_rt`, short-XIC estimator); null on edge/short traces |
| `ms1_rt_fit_status` | str | `ok` / failure reason / `not_extracted` |
| `ms1_im_apex/fwhm/sigma/skew/r2` | float64 | IM peak fit (1/K0) |
| `ms1_im_snr` | float64 | IM mobilogram SNR (`trace_snr_v1`) |
| `ms1_im_fit_status` | str | as above |
| `rt_width_reliable` / `im_width_reliable` / `width_reliable` | bool | ⚠️ **provisional** band (`r2≥0.8 & snr≥20`); gate width labels on these |

### Isotope, intensity, counts
| field | type | notes |
|---|---|---|
| `isotope_cosim` | float64 | isotope-envelope cosine similarity |
| `ms1_iso_0`…`ms1_iso_4` | float64 | isotope-peak intensities |
| `precursor_intensity`, `ms1_total_intensity`, `fragment_total_intensity` | float64 | |
| `n_peaks`, `n_fragments_merged` | int64 | |

### Derived presets
| field | type | notes |
|---|---|---|
| `strict` | bool | `n_engines==2 & both q≤0.01 & width_reliable` — one-click high-confidence subset |
| `split` | str | `train` / `validation` / `test` |

---

## `tier3_fragments` — one row per matched b/y fragment

Grain: `(raw_file, precursor_id, fragment_type, fragment_ordinal, fragment_charge)`
of the canonical peptidoform. Fragments are imspy-rematched for both engines
under identical settings (deterministic one-to-one matching, 20 ppm).

| field | type | notes |
|---|---|---|
| `accession`, `raw_file`, `precursor_id`, `sage_psm_id` | | join keys back to Tier 1 |
| `assignment` | str | `consensus` (both engines agree) / `sage` / `fragpipe` |
| `canonical_engine` | str | engine whose peptide was matched (`sage` preferred) |
| `sequence`, `modified_sequence`, `sequence_normalized` | str | the matched peptidoform |
| `precursor_charge` | int64 | |
| `n_engines`, `peptidoform_conflict` | int64 / bool | engine agreement context |
| `fragment_type` | str | `b` / `y` |
| `fragment_ordinal` | int64 | 1…N−1 |
| `fragment_charge` | int64 | fragment charge |
| `fragment_mz_calculated` | float64 | theoretical m/z |
| `fragment_mz_experimental` | float64 | matched peak m/z |
| `fragment_intensity` | float64 | matched peak intensity (the prediction target) |
| `ppm_error` | float64 | `(exp − calc)/calc · 1e6` |
| `match_method` | str | `imspy_rematch` |
| `sage_native_matched` | bool | also reported by Sage's own annotator (provenance) |
| `sage_native_intensity` | float64 | Sage-native intensity (null if not natively matched) |
| `split` | str | `train` / `validation` / `test` |

---

*Schema versions:* tier1 `tier1.v5`, tier3 `tier3.v1`. See `manifest.json` for the
build-pipeline commit, filter parameters, and per-dataset provenance.
