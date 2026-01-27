#!/usr/bin/env python3
"""
Train CCS prediction model using imspy-predictors.

Uses the ionmob-style architecture (GRU-based) for CCS prediction from peptide sequences.
"""

import json
import sys
from pathlib import Path
from typing import Tuple, Dict, Any

import numpy as np
import pandas as pd
import yaml

# imspy-predictors imports
from imspy_predictors.ccs import load_deep_ccs_predictor, DeepCCSPredictor
from imspy_predictors.utilities.tokenizer import PeptideTokenizer

# TensorFlow imports
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
from tensorflow.keras.callbacks import EarlyStopping, ModelCheckpoint, ReduceLROnPlateau

# Snakemake provides these variables
data_file = Path(snakemake.input.data)
model_output = Path(snakemake.output.model)
metrics_output = Path(snakemake.output.metrics)
config_output = Path(snakemake.output.config)
log_file = Path(snakemake.log[0])

# Training parameters
test_split = snakemake.params.test_split
val_split = snakemake.params.val_split
batch_size = snakemake.params.batch_size
epochs = snakemake.params.epochs
patience = snakemake.params.patience


def setup_logging(log_path: Path):
    """Set up logging to file and stdout."""
    import logging

    log_path.parent.mkdir(parents=True, exist_ok=True)

    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)

    fh = logging.FileHandler(log_path)
    fh.setLevel(logging.INFO)

    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)

    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)

    logger.addHandler(fh)
    logger.addHandler(ch)

    return logger


def load_training_data(data_path: Path) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Load and prepare training data.

    Returns train, validation, and test DataFrames.
    """
    logger.info(f"Loading data from {data_path}")
    df = pd.read_parquet(data_path)

    # Require necessary columns
    required_cols = ['sequence', 'charge', 'ccs']
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    # Filter valid entries
    df = df[df['ccs'].notna() & (df['ccs'] > 0)]
    df = df[df['sequence'].notna() & (df['sequence'] != '')]
    df = df[df['charge'].notna() & (df['charge'] > 0)]

    logger.info(f"Valid entries after filtering: {len(df)}")

    # Use pre-computed split if available
    if 'split' in df.columns:
        train_df = df[df['split'] == 'train'].copy()
        test_df = df[df['split'] == 'test'].copy()
    else:
        # Create split by peptide
        unique_peptides = df['sequence'].unique()
        n_test = int(len(unique_peptides) * test_split)

        np.random.seed(42)
        test_peptides = set(np.random.choice(unique_peptides, n_test, replace=False))

        test_df = df[df['sequence'].isin(test_peptides)].copy()
        train_df = df[~df['sequence'].isin(test_peptides)].copy()

    # Further split train into train/validation
    train_peptides = train_df['sequence'].unique()
    n_val = int(len(train_peptides) * val_split)

    np.random.seed(42)
    val_peptides = set(np.random.choice(train_peptides, n_val, replace=False))

    val_df = train_df[train_df['sequence'].isin(val_peptides)].copy()
    train_df = train_df[~train_df['sequence'].isin(val_peptides)].copy()

    logger.info(f"Train: {len(train_df)} records ({train_df['sequence'].nunique()} peptides)")
    logger.info(f"Validation: {len(val_df)} records ({val_df['sequence'].nunique()} peptides)")
    logger.info(f"Test: {len(test_df)} records ({test_df['sequence'].nunique()} peptides)")

    return train_df, val_df, test_df


def create_tokenizer() -> PeptideTokenizer:
    """Create peptide tokenizer for sequence encoding."""
    # Standard amino acids + common modifications
    tokenizer = PeptideTokenizer()
    return tokenizer


def prepare_data_for_training(
    df: pd.DataFrame,
    tokenizer: PeptideTokenizer,
    max_length: int = 50,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Prepare data for model training.

    Returns:
        sequences: Encoded peptide sequences (N, max_length)
        charges: Charge states (N,)
        ccs_values: Target CCS values (N,)
    """
    sequences = df['sequence'].values
    charges = df['charge'].values.astype(np.int32)
    ccs_values = df['ccs'].values.astype(np.float32)

    # Encode sequences using tokenizer
    encoded_sequences = []
    for seq in sequences:
        encoded = tokenizer.encode(seq, max_length=max_length)
        encoded_sequences.append(encoded)

    sequences_array = np.array(encoded_sequences, dtype=np.int32)

    return sequences_array, charges, ccs_values


def build_ccs_model(
    vocab_size: int = 30,
    max_length: int = 50,
    embedding_dim: int = 128,
    gru_units: int = 256,
    dropout_rate: float = 0.2,
) -> keras.Model:
    """
    Build the CCS prediction model.

    Architecture based on ionmob:
    - Embedding layer for amino acid sequences
    - Bidirectional GRU layers
    - Dense layers with charge state input
    - Linear regression head for CCS prediction
    """
    # Sequence input
    sequence_input = keras.Input(shape=(max_length,), dtype='int32', name='sequence')

    # Charge input
    charge_input = keras.Input(shape=(1,), dtype='float32', name='charge')

    # Embedding layer
    x = layers.Embedding(vocab_size, embedding_dim, mask_zero=True)(sequence_input)

    # Bidirectional GRU layers
    x = layers.Bidirectional(layers.GRU(gru_units, return_sequences=True, dropout=dropout_rate))(x)
    x = layers.Bidirectional(layers.GRU(gru_units // 2, dropout=dropout_rate))(x)

    # Concatenate with charge
    x = layers.Concatenate()([x, charge_input])

    # Dense layers
    x = layers.Dense(256, activation='relu')(x)
    x = layers.Dropout(dropout_rate)(x)
    x = layers.Dense(128, activation='relu')(x)
    x = layers.Dropout(dropout_rate)(x)
    x = layers.Dense(64, activation='relu')(x)

    # Output layer (CCS prediction)
    output = layers.Dense(1, activation='linear', name='ccs')(x)

    model = keras.Model(inputs=[sequence_input, charge_input], outputs=output)

    return model


def train_model(
    model: keras.Model,
    train_data: Tuple[np.ndarray, np.ndarray, np.ndarray],
    val_data: Tuple[np.ndarray, np.ndarray, np.ndarray],
    model_path: Path,
    batch_size: int,
    epochs: int,
    patience: int,
) -> Dict[str, Any]:
    """Train the CCS model."""
    logger.info("Starting model training...")

    X_train_seq, X_train_charge, y_train = train_data
    X_val_seq, X_val_charge, y_val = val_data

    # Compile model
    model.compile(
        optimizer=keras.optimizers.Adam(learning_rate=1e-3),
        loss='mse',
        metrics=['mae']
    )

    logger.info(f"Model parameters: {model.count_params():,}")
    model.summary(print_fn=logger.info)

    # Callbacks
    callbacks = [
        EarlyStopping(
            monitor='val_loss',
            patience=patience,
            restore_best_weights=True,
            verbose=1
        ),
        ModelCheckpoint(
            str(model_path),
            monitor='val_loss',
            save_best_only=True,
            verbose=1
        ),
        ReduceLROnPlateau(
            monitor='val_loss',
            factor=0.5,
            patience=patience // 2,
            min_lr=1e-6,
            verbose=1
        ),
    ]

    # Train
    history = model.fit(
        x={'sequence': X_train_seq, 'charge': X_train_charge.reshape(-1, 1)},
        y=y_train,
        validation_data=(
            {'sequence': X_val_seq, 'charge': X_val_charge.reshape(-1, 1)},
            y_val
        ),
        batch_size=batch_size,
        epochs=epochs,
        callbacks=callbacks,
        verbose=1
    )

    return history.history


def evaluate_model(
    model: keras.Model,
    test_data: Tuple[np.ndarray, np.ndarray, np.ndarray],
) -> Dict[str, float]:
    """Evaluate model on test set."""
    logger.info("Evaluating model on test set...")

    X_test_seq, X_test_charge, y_test = test_data

    # Predict
    y_pred = model.predict(
        {'sequence': X_test_seq, 'charge': X_test_charge.reshape(-1, 1)},
        verbose=0
    ).flatten()

    # Calculate metrics
    from sklearn.metrics import r2_score, mean_absolute_error, mean_squared_error

    r2 = r2_score(y_test, y_pred)
    mae = mean_absolute_error(y_test, y_pred)
    rmse = np.sqrt(mean_squared_error(y_test, y_pred))

    # Relative errors
    relative_errors = np.abs((y_test - y_pred) / y_test) * 100
    mape = np.mean(relative_errors)
    median_relative_error = np.median(relative_errors)

    # Percentage within thresholds
    within_1pct = (relative_errors < 1.0).mean() * 100
    within_2pct = (relative_errors < 2.0).mean() * 100
    within_5pct = (relative_errors < 5.0).mean() * 100

    metrics = {
        'r2': float(r2),
        'mae': float(mae),
        'rmse': float(rmse),
        'mape': float(mape),
        'median_relative_error': float(median_relative_error),
        'within_1pct': float(within_1pct),
        'within_2pct': float(within_2pct),
        'within_5pct': float(within_5pct),
        'n_test': len(y_test),
        'ccs_mean': float(y_test.mean()),
        'ccs_std': float(y_test.std()),
        'ccs_min': float(y_test.min()),
        'ccs_max': float(y_test.max()),
    }

    logger.info("=" * 50)
    logger.info("Test Set Evaluation Results")
    logger.info("=" * 50)
    logger.info(f"R²: {r2:.4f}")
    logger.info(f"MAE: {mae:.2f} Å²")
    logger.info(f"RMSE: {rmse:.2f} Å²")
    logger.info(f"MAPE: {mape:.2f}%")
    logger.info(f"Median Relative Error: {median_relative_error:.2f}%")
    logger.info(f"Within 1%: {within_1pct:.1f}%")
    logger.info(f"Within 2%: {within_2pct:.1f}%")
    logger.info(f"Within 5%: {within_5pct:.1f}%")
    logger.info("=" * 50)
    logger.info(f"Target (Meier et al.): R > 0.99, median error < 1.4%")
    logger.info("=" * 50)

    return metrics


def main():
    global logger
    logger = setup_logging(log_file)

    logger.info("=" * 60)
    logger.info("CLAUDIUS-PROTEOMICS: CCS Model Training")
    logger.info("=" * 60)

    # Create output directories
    model_output.parent.mkdir(parents=True, exist_ok=True)

    # Set random seeds for reproducibility
    np.random.seed(42)
    tf.random.set_seed(42)

    # Load data
    train_df, val_df, test_df = load_training_data(data_file)

    # Create tokenizer
    tokenizer = create_tokenizer()
    vocab_size = tokenizer.vocab_size

    logger.info(f"Tokenizer vocab size: {vocab_size}")

    # Prepare data
    logger.info("Preparing training data...")
    train_data = prepare_data_for_training(train_df, tokenizer)
    val_data = prepare_data_for_training(val_df, tokenizer)
    test_data = prepare_data_for_training(test_df, tokenizer)

    logger.info(f"Train sequences shape: {train_data[0].shape}")
    logger.info(f"Val sequences shape: {val_data[0].shape}")
    logger.info(f"Test sequences shape: {test_data[0].shape}")

    # Build model
    logger.info("Building CCS prediction model...")
    model = build_ccs_model(
        vocab_size=vocab_size,
        max_length=50,
        embedding_dim=128,
        gru_units=256,
        dropout_rate=0.2,
    )

    # Train
    history = train_model(
        model,
        train_data,
        val_data,
        model_output,
        batch_size,
        epochs,
        patience,
    )

    # Load best model
    logger.info(f"Loading best model from {model_output}")
    model = keras.models.load_model(model_output)

    # Evaluate
    metrics = evaluate_model(model, test_data)
    metrics['training_history'] = {
        'loss': [float(x) for x in history.get('loss', [])],
        'val_loss': [float(x) for x in history.get('val_loss', [])],
        'mae': [float(x) for x in history.get('mae', [])],
        'val_mae': [float(x) for x in history.get('val_mae', [])],
    }

    # Save metrics
    logger.info(f"Saving metrics to {metrics_output}")
    with open(metrics_output, 'w') as f:
        json.dump(metrics, f, indent=2)

    # Save training config
    config = {
        'test_split': test_split,
        'val_split': val_split,
        'batch_size': batch_size,
        'epochs': epochs,
        'patience': patience,
        'vocab_size': vocab_size,
        'max_length': 50,
        'embedding_dim': 128,
        'gru_units': 256,
        'dropout_rate': 0.2,
        'n_train': len(train_df),
        'n_val': len(val_df),
        'n_test': len(test_df),
        'n_train_peptides': int(train_df['sequence'].nunique()),
        'n_val_peptides': int(val_df['sequence'].nunique()),
        'n_test_peptides': int(test_df['sequence'].nunique()),
    }

    logger.info(f"Saving config to {config_output}")
    with open(config_output, 'w') as f:
        yaml.dump(config, f, default_flow_style=False)

    logger.info("=" * 60)
    logger.info("Training completed successfully!")
    logger.info(f"Model saved to: {model_output}")
    logger.info(f"Metrics saved to: {metrics_output}")
    logger.info("=" * 60)


if __name__ == "__main__":
    main()
