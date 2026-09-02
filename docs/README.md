# Documentation index

## Start here
| Document | What it covers |
|---|---|
| [DESIGN.md](DESIGN.md) | The full project design — vision, corpus principles, architecture. |
| [SAN_JOSE_PITCH.md](SAN_JOSE_PITCH.md) | Short-form motivation for the "San José" reference layer. |
| [RUNNER_ARCHITECTURE.md](RUNNER_ARCHITECTURE.md) | How the 6-step per-dataset pipeline runner is put together. |

## Corpus & data model
| Document | What it covers |
|---|---|
| [CORPUS_SCHEMA.md](CORPUS_SCHEMA.md) | Field-by-field reference for the published parquet corpus. |
| [DATASET_CARD.md](DATASET_CARD.md) | Hugging Face dataset card for the released corpus. |
| [DATASET_DEFINITION.md](DATASET_DEFINITION.md) | What counts as one dataset, and how runs are grouped. |
| [RUNNER_OUTPUT_SCHEMA.md](RUNNER_OUTPUT_SCHEMA.md) | Per-step output contracts of the runner. |
| [PRECURSOR_BLOB_DESIGN.md](PRECURSOR_BLOB_DESIGN.md) | Layout of the per-precursor raw-signal blobs. |
| [RAW_DATA_ARCHITECTURE.md](RAW_DATA_ARCHITECTURE.md) | How Bruker `.d` data is staged, normalised and addressed. |

## Selection, provenance & validation
| Document | What it covers |
|---|---|
| [DATASET_PRIORITY_LIST.md](DATASET_PRIORITY_LIST.md) | Which PRIDE datasets were selected, and on what criteria. |
| [RUN_LEDGER.md](RUN_LEDGER.md) | Provenance: code, config, job and engine versions behind each result. |
| [HF_CORPUS_VALIDATION.md](HF_CORPUS_VALIDATION.md) | End-to-end model training used to prove the corpus carries signal. |
| [PRE_SCALE_DECISIONS.md](PRE_SCALE_DECISIONS.md) | Decisions locked in before scaling processing out. |

## Modelling notes
| Document | What it covers |
|---|---|
| [LAB_DRIFT_LATENT_MODEL.md](LAB_DRIFT_LATENT_MODEL.md) | Latent-variable treatment of lab/batch drift. |

## Figures
`pipeline_dag.{dot,svg,png}` — rendered Snakemake DAG of the processing pipeline.
