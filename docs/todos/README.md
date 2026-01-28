# San José Database - TODO Distribution

This directory contains work packages for the San José peptide property database project.

## Work Packages

| # | File | Area | Priority | Est. Complexity |
|---|------|------|----------|-----------------|
| 01 | [raw-data-extraction.md](01-raw-data-extraction.md) | RAW DATA Extraction | Medium | Medium |
| 02 | [progress-monitoring.md](02-progress-monitoring.md) | Progress Monitoring | Medium | Medium |
| 03 | [san-jose-database-design.md](03-san-jose-database-design.md) | Database Design | Critical | High |
| 04 | [third-party-id-fields.md](04-third-party-id-fields.md) | ID Field Extraction | Critical | Medium |
| 05 | [overlap-union-datasets.md](05-overlap-union-datasets.md) | Consensus Strategy | High | Medium |
| 06 | [model-grid-search.md](06-model-grid-search.md) | Hyperparameter Optimization | Low | High |

## Dependencies

```
03 (Database Design) ──► 04 (ID Fields) ──► 05 (Overlap/Union)
                                                    │
                      ┌─────────────────────────────┘
                      ▼
              01 (RAW Extraction)
                      │
                      ▼
              02 (Monitoring)
                      │
                      ▼
              06 (Grid Search)
```

**Recommended order:** 03 → 04 → 05 → 01 → 02 → 06

**Rationale:**
- 04 (ID Fields) must come before 05 - consensus matching needs q-values and PEP scores
- 01 (RAW Extraction) can wait - basic CCS/RT from search engines is sufficient initially
- 06 (Grid Search) is last - need substantial data before hyperparameter optimization matters

## Assignment Tracking

| Package | Owner | Status | Start Date | Target Date |
|---------|-------|--------|------------|-------------|
| 01 | TBD | Not Started | | |
| 02 | TBD | Not Started | | |
| 03 | TBD | Not Started | | |
| 04 | TBD | Not Started | | |
| 05 | TBD | Not Started | | |
| 06 | TBD | Not Started | | |

## Context

- **Project:** CLAUDIUS-PROTEOMICS
- **Database Name:** San José
- **Goal:** Build the world's largest peptide property prediction resource
- **Data Source:** PRIDE Archive (public proteomics data)
- **Search Engines:** FragPipe + DIA-NN (orthogonal validation)
