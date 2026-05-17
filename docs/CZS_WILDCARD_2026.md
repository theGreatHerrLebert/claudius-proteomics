# CZS Wildcard 2026 — 1-Page Proposal

---

## The Scientist Who Shouldn't Exist

**One PhD student. One AI. An infrastructure that should have taken a team of ten.**

---

### The Observation

Over the past 12 months, a single doctoral researcher built *San José* — a fully operational pipeline that systematically reprocesses petabytes of public mass spectrometry data from the PRIDE repository through three independent search engines, extracts raw 4D signal features, and produces bias-aware, versioned training datasets for machine learning models in proteomics.

This pipeline includes: automated data download, triple-engine search (FragPipe, DIA-NN, Sage), consensus stratification, raw signal extraction from Bruker timsTOF instruments, an interactive 4D data browser (React + FastAPI + Deck.gl), collision energy calibration, CCS reproducibility analysis, and a distributable archive format. The codebase spans ~15,000 lines across Python, Rust, TypeScript, and Snakemake.

**This should not have been possible.** Not by one person, not in one year.

It was built entirely through *conversational programming with AI* — a practice the developer community calls **"vibe coding"**.

---

### The Wild Question

> **What happens to science when the engineering cost of building research infrastructure drops to near zero?**

We propose to find out — rigorously.

Today, the bottleneck in computational science is rarely the algorithm. It is the **engineering**: the data pipelines, format conversions, database schemas, visualizations, deployment scripts, and test harnesses that make research *reproducible and reusable*. This infrastructure work is invisible, unglamorous, and traditionally requires dedicated software engineering teams that most labs cannot afford.

AI-assisted coding fundamentally changes this equation. But no one has systematically studied *how*, *when it works*, *when it fails*, and *what it means for how we organize scientific research*.

---

### What We Propose

**A 2-year, interdisciplinary study of AI-assisted scientific software development ("vibe coding") as a new research methodology**, using proteomics infrastructure as a living laboratory.

**Three workstreams:**

1. **Empirical methodology study** — We will document, measure, and analyze the AI-assisted development process across multiple scientific domains. What types of tasks benefit most? Where does AI-assisted coding introduce subtle errors? How does code quality compare to traditionally developed scientific software? We will develop a taxonomy of failure modes and a practical framework for "when to trust the AI and when to verify."

2. **Stress-test at scale** — We will push *San José* from proof-of-concept to production: reprocessing 100+ PRIDE datasets (~500 TB raw data), validating against published benchmarks (Meier et al. 2021: CCS prediction R > 0.99), and training open models for CCS, retention time, and fragmentation prediction. This serves as both a scientific deliverable and a real-world testbed for the methodology.

3. **Transferability experiment** — Can the approach transfer? We will replicate the "one researcher + AI builds infrastructure" model in [2–3 additional scientific domains, e.g., cryo-EM data processing, genomic variant calling, climate model post-processing]. If the approach generalizes, the implications are profound. If it doesn't, understanding *why not* is equally valuable.

---

### Why This Is Wild

- **It challenges the team model.** Modern science assumes that complex software requires dedicated engineering teams. We claim this assumption is becoming obsolete — and we have the proof-of-concept to back it up.
- **It is self-referential.** We propose to use AI-assisted coding to study AI-assisted coding. The methodology is its own subject.
- **It is uncomfortable.** If one person + AI can build what previously required ten, this has implications for hiring, funding structures, and the division of labor in computational science. We do not shy away from these questions.
- **It could be wrong.** Perhaps San José succeeded despite vibe coding, not because of it. Perhaps the approach only works for certain types of infrastructure. The honest answer is: we don't know yet. That's why it needs rigorous study.

---

### The Team

| Role | Expertise | Contribution |
|------|-----------|-------------|
| **[PI 1]** | Proteomics / Mass Spectrometry | Domain science, San José pipeline, vibe coding practitioner |
| **[PI 2]** | Software Engineering / Empirical SE | Methodology design, code quality analysis, failure taxonomy |
| **[PI 3]** | [Third discipline, e.g., Science of Science / Research Policy] | Transferability study, implications for research organization |

---

### Budget Sketch (€750k / 2 years)

| Item | Amount | Purpose |
|------|--------|---------|
| 2 Postdocs (E13, 24 months) | ~€400k | Methodology study + transferability experiments |
| 1 PhD student top-up | ~€50k | San José scale-up and validation |
| Cloud compute (HPC, GPU) | ~€150k | Large-scale reprocessing of PRIDE data |
| Travel, workshops, open access | ~€50k | Dissemination, community building |
| Equipment / licenses | ~€100k | AI coding tool subscriptions, storage |

*+20% overhead (CZS Pauschale)*

---

### Expected Outcomes

1. **A practical framework** for AI-assisted scientific software development — when it works, when it fails, how to verify
2. **San José v1.0** — the first open, bias-aware, multi-engine reference layer for ion mobility proteomics (public dataset + trained models)
3. **Transferability evidence** — does the "one researcher + AI" model generalize beyond proteomics?
4. **An honest assessment** of what this means for how we fund, organize, and evaluate computational research

---

*"The best way to predict the future is to build it — and then study what you built."*
