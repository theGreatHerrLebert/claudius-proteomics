# Model Grid Search Strategy

**Owner:** TBD
**Priority:** Medium
**Status:** Not Started

## Overview

Design and implement a systematic hyperparameter optimization strategy for deep learning models (CCS prediction, RT prediction, etc.). Use modern optimization techniques to find optimal architectures and training configurations.

## Current Model

- **Architecture:** ionmob-style GRU for CCS prediction
- **Framework:** TensorFlow 2.16+
- **Location:** `scripts/train_ccs.py`

## Hyperparameters to Optimize

### 1. Architecture Parameters

| Parameter | Search Space | Default |
|-----------|--------------|---------|
| RNN type | GRU, LSTM, Bidirectional-GRU, Bidirectional-LSTM | GRU |
| Number of layers | 1, 2, 3 | 2 |
| Hidden units | 64, 128, 256, 512 | 128 |
| Embedding dimension | 16, 32, 64, 128 | 64 |
| Dropout rate | 0.0, 0.1, 0.2, 0.3, 0.5 | 0.2 |
| Attention mechanism | None, Self-Attention, Multi-Head | None |
| Dense layers after RNN | 0, 1, 2 | 1 |
| Dense layer units | 32, 64, 128 | 64 |

### 2. Training Parameters

| Parameter | Search Space | Default |
|-----------|--------------|---------|
| Learning rate | 1e-4, 5e-4, 1e-3, 5e-3 | 1e-3 |
| Batch size | 64, 128, 256, 512, 1024 | 256 |
| Optimizer | Adam, AdamW, SGD+momentum | Adam |
| LR schedule | Constant, Cosine, ReduceOnPlateau | ReduceOnPlateau |
| Weight decay | 0, 1e-5, 1e-4, 1e-3 | 0 |
| Gradient clipping | None, 1.0, 5.0 | None |
| Epochs | 50, 100, 200 | 100 |

### 3. Data/Input Parameters

| Parameter | Search Space | Default |
|-----------|--------------|---------|
| Sequence encoding | One-hot, Learned embedding, Pre-trained (ESM-2) | Learned |
| Max sequence length | 30, 50, 100 | 50 |
| Charge encoding | Ordinal, One-hot, Embedding | One-hot |
| Include mass feature | Yes, No | No |
| Include RT feature | Yes, No | No |
| Data augmentation | None, Reverse, AA-swap | None |

### 4. Regularization Parameters

| Parameter | Search Space | Default |
|-----------|--------------|---------|
| Early stopping patience | 5, 10, 20 | 10 |
| Label smoothing | 0, 0.05, 0.1 | 0 |
| Batch normalization | Yes, No | No |
| Layer normalization | Yes, No | No |

## Search Strategy

### Phase 1: Random Search (Broad Exploration)

**Goal:** Identify promising regions of hyperparameter space.

```python
n_trials = 100
epochs_per_trial = 20
strategy = "random"
```

- Sample 100 random configurations
- Train each for 20 epochs
- Use validation MAE as objective
- Identify top 10% configurations

### Phase 2: Bayesian Optimization (Fine-Tuning)

**Goal:** Efficiently explore promising regions.

```python
n_trials = 200
strategy = "TPE"  # Tree-structured Parzen Estimator
pruner = "Hyperband"  # Early stopping of bad trials
```

- Use Optuna with TPE sampler
- Hyperband pruning (stop bad trials early)
- Focus on promising regions from Phase 1

### Phase 3: Final Training

**Goal:** Validate top configurations with full training.

```python
top_k = 3
epochs = 200
cross_validation_folds = 5
```

- Take top 3 configurations
- Full training (200 epochs, early stopping)
- 5-fold cross-validation for variance estimation
- Report mean ± std for all metrics

## Implementation

### Optuna Integration

```python
# scripts/hyperparameter_search.py

import optuna
from optuna.integration import TFKerasPruningCallback

def objective(trial):
    # Sample hyperparameters
    params = {
        'rnn_type': trial.suggest_categorical('rnn_type', ['GRU', 'LSTM']),
        'n_layers': trial.suggest_int('n_layers', 1, 3),
        'hidden_units': trial.suggest_categorical('hidden_units', [64, 128, 256]),
        'embedding_dim': trial.suggest_categorical('embedding_dim', [32, 64, 128]),
        'dropout': trial.suggest_float('dropout', 0.0, 0.5),
        'learning_rate': trial.suggest_float('learning_rate', 1e-4, 1e-2, log=True),
        'batch_size': trial.suggest_categorical('batch_size', [128, 256, 512]),
    }

    # Build model
    model = build_model(**params)

    # Train with pruning callback
    history = model.fit(
        train_data,
        validation_data=val_data,
        epochs=50,
        callbacks=[
            TFKerasPruningCallback(trial, 'val_mae'),
            EarlyStopping(patience=5)
        ]
    )

    return min(history.history['val_mae'])

# Run optimization
study = optuna.create_study(
    direction='minimize',
    pruner=optuna.pruners.HyperbandPruner(),
    storage='sqlite:///models/grid_search/study.db'
)
study.optimize(objective, n_trials=200)
```

### Snakemake Rule

```python
# rules/grid_search.smk

rule hyperparameter_search:
    input:
        snapshot="database/snapshots/{version}/train.parquet"
    output:
        study="models/grid_search/{version}/study.db",
        best_params="models/grid_search/{version}/best_params.yaml"
    params:
        n_trials=config.get("grid_search_trials", 200)
    resources:
        gpu=1,
        time="48:00:00"
    script:
        "../scripts/hyperparameter_search.py"

rule analyze_grid_search:
    input:
        study="models/grid_search/{version}/study.db"
    output:
        importance="models/grid_search/{version}/importance.png",
        history="models/grid_search/{version}/optimization.png",
        report="models/grid_search/{version}/report.html"
    script:
        "../scripts/analyze_grid_search.py"
```

## Resource Requirements

### Compute Budget

| Phase | Trials | Epochs/Trial | GPU Hours (est.) |
|-------|--------|--------------|------------------|
| Random Search | 100 | 20 | ~50 |
| Bayesian Opt | 200 | 50 (avg with pruning) | ~200 |
| Final Training | 3×5 | 200 | ~75 |
| **Total** | | | **~325** |

### HPC Configuration

```yaml
# profiles/mogon2/config.yaml - grid_search resources
set-resources:
  hyperparameter_search:
    slurm_partition: "gpu"
    slurm_extra: "'--gres=gpu:1'"
    mem_mb: 32000
    time: "48:00:00"
```

## Output Artifacts

```
models/grid_search/{version}/
├── study.db              # Optuna study database (SQLite)
├── best_params.yaml      # Optimal hyperparameters
├── trials.parquet        # All trial results
├── importance.png        # Hyperparameter importance plot
├── optimization.png      # Optimization history plot
├── parallel_coord.png    # Parallel coordinates plot
├── slice.png             # Slice plot for top params
└── report.html           # Comprehensive HTML report
```

## Analysis Outputs

### Hyperparameter Importance

Identify which parameters matter most:

```python
importance = optuna.importance.get_param_importances(study)
# Expected: learning_rate > hidden_units > dropout > ...
```

### Visualization

1. **Optimization History:** MAE vs trial number
2. **Parallel Coordinates:** Parameter combinations for top trials
3. **Slice Plot:** Individual parameter effects
4. **Contour Plot:** Interaction between top 2 parameters

## Alternative: Ray Tune

For distributed search across multiple GPUs:

```python
from ray import tune
from ray.tune.schedulers import ASHAScheduler

scheduler = ASHAScheduler(
    max_t=100,
    grace_period=10,
    reduction_factor=3
)

analysis = tune.run(
    train_function,
    config=search_space,
    num_samples=200,
    scheduler=scheduler,
    resources_per_trial={"cpu": 4, "gpu": 1}
)
```

## Acceptance Criteria

- [ ] Optuna integrated with TensorFlow training
- [ ] Hyperband pruning working (bad trials stopped early)
- [ ] Study persisted to SQLite for resume capability
- [ ] Best parameters exported to YAML
- [ ] Importance and history plots generated
- [ ] Top 3 configs validated with cross-validation
- [ ] Final model trained with optimal hyperparameters
- [ ] Documentation of search results and recommendations
