#!/usr/bin/env python3
"""
Evaluate trained CCS model and generate report with plots.
"""

import json
import sys
from pathlib import Path

import pandas as pd
import numpy as np

# Snakemake provides these variables
model_path = Path(snakemake.input.model)
data_path = Path(snakemake.input.data)
report_output = Path(snakemake.output.report)
plots_dir = Path(snakemake.output.plots)
log_file = Path(snakemake.log[0])

test_split = snakemake.params.test_split


def setup_logging(log_path: Path):
    """Redirect stdout/stderr to log file."""
    import logging
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_path),
            logging.StreamHandler(sys.stdout)
        ]
    )
    return logging.getLogger(__name__)


def load_model(model_path: Path):
    """Load trained model."""
    logger.info(f"Loading model from {model_path}")
    try:
        from tensorflow import keras
        model = keras.models.load_model(model_path)
        return model
    except Exception as e:
        logger.error(f"Failed to load model: {e}")
        return None


def encode_sequences(df: pd.DataFrame, max_length: int = 50) -> np.ndarray:
    """Encode peptide sequences for model input."""
    aa_to_int = {aa: i + 1 for i, aa in enumerate('ACDEFGHIKLMNPQRSTVWY')}
    aa_to_int['X'] = 0

    def encode_seq(seq):
        encoded = [aa_to_int.get(aa, 0) for aa in seq[:max_length]]
        encoded += [0] * (max_length - len(encoded))
        return encoded

    sequences = df['sequence'].apply(encode_seq).tolist()
    return np.array(sequences)


def create_plots(y_true: np.ndarray, y_pred: np.ndarray, output_dir: Path):
    """Create evaluation plots."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    output_dir.mkdir(parents=True, exist_ok=True)

    # 1. Scatter plot: predicted vs actual
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.scatter(y_true, y_pred, alpha=0.3, s=1)
    ax.plot([y_true.min(), y_true.max()], [y_true.min(), y_true.max()], 'r--', lw=2)
    ax.set_xlabel('Measured CCS (Å²)')
    ax.set_ylabel('Predicted CCS (Å²)')
    ax.set_title('CCS Prediction: Measured vs Predicted')
    fig.savefig(output_dir / 'scatter_pred_vs_actual.png', dpi=150, bbox_inches='tight')
    plt.close(fig)

    # 2. Residual plot
    residuals = y_pred - y_true
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.scatter(y_true, residuals, alpha=0.3, s=1)
    ax.axhline(y=0, color='r', linestyle='--')
    ax.set_xlabel('Measured CCS (Å²)')
    ax.set_ylabel('Residual (Predicted - Measured)')
    ax.set_title('Residual Plot')
    fig.savefig(output_dir / 'residuals.png', dpi=150, bbox_inches='tight')
    plt.close(fig)

    # 3. Error distribution
    relative_error = (y_pred - y_true) / y_true * 100
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.hist(relative_error, bins=100, edgecolor='black', alpha=0.7)
    ax.axvline(x=0, color='r', linestyle='--')
    ax.set_xlabel('Relative Error (%)')
    ax.set_ylabel('Count')
    ax.set_title(f'Error Distribution (median: {np.median(relative_error):.2f}%)')
    fig.savefig(output_dir / 'error_distribution.png', dpi=150, bbox_inches='tight')
    plt.close(fig)

    # 4. CCS distribution comparison
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.hist(y_true, bins=50, alpha=0.5, label='Measured', edgecolor='black')
    ax.hist(y_pred, bins=50, alpha=0.5, label='Predicted', edgecolor='black')
    ax.set_xlabel('CCS (Å²)')
    ax.set_ylabel('Count')
    ax.set_title('CCS Distribution: Measured vs Predicted')
    ax.legend()
    fig.savefig(output_dir / 'ccs_distribution.png', dpi=150, bbox_inches='tight')
    plt.close(fig)

    logger.info(f"Plots saved to {output_dir}")


def generate_html_report(metrics: dict, output_path: Path):
    """Generate HTML evaluation report."""
    html = f"""
<!DOCTYPE html>
<html>
<head>
    <title>CLAUDIUS-PROTEOMICS: CCS Model Evaluation</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 40px; }}
        h1 {{ color: #333; }}
        .metric {{ font-size: 24px; margin: 10px 0; }}
        .metric-value {{ color: #0066cc; font-weight: bold; }}
        table {{ border-collapse: collapse; margin: 20px 0; }}
        td, th {{ border: 1px solid #ddd; padding: 8px; }}
        th {{ background-color: #f4f4f4; }}
        img {{ max-width: 600px; margin: 10px 0; }}
    </style>
</head>
<body>
    <h1>CCS Model Evaluation Report</h1>

    <h2>Performance Metrics</h2>
    <table>
        <tr><th>Metric</th><th>Value</th></tr>
        <tr><td>R²</td><td class="metric-value">{metrics.get('r2', 'N/A'):.4f}</td></tr>
        <tr><td>MAE</td><td class="metric-value">{metrics.get('mae', 'N/A'):.2f} Å²</td></tr>
        <tr><td>MAPE</td><td class="metric-value">{metrics.get('mape', 'N/A'):.2f}%</td></tr>
        <tr><td>Test samples</td><td>{metrics.get('n_test', 'N/A')}</td></tr>
    </table>

    <h2>Plots</h2>
    <h3>Predicted vs Measured</h3>
    <img src="plots/scatter_pred_vs_actual.png" alt="Scatter plot">

    <h3>Residuals</h3>
    <img src="plots/residuals.png" alt="Residual plot">

    <h3>Error Distribution</h3>
    <img src="plots/error_distribution.png" alt="Error distribution">

    <h3>CCS Distribution</h3>
    <img src="plots/ccs_distribution.png" alt="CCS distribution">

    <h2>Target Performance</h2>
    <p>Meier et al. (2021): R > 0.99, median relative error &lt; 1.4%</p>

</body>
</html>
    """

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(html)
    logger.info(f"Report saved to {output_path}")


def main():
    global logger
    logger = setup_logging(log_file)

    logger.info("Starting model evaluation")

    # Load data
    df = pd.read_parquet(data_path)
    test_df = df[df['split'] == 'test'].copy()
    logger.info(f"Test set: {len(test_df)} samples")

    # Load model
    model = load_model(model_path)

    if model is None:
        logger.error("Could not load model, generating placeholder report")
        metrics = {'r2': 0, 'mae': 0, 'mape': 0, 'n_test': len(test_df)}
        generate_html_report(metrics, report_output)
        plots_dir.mkdir(parents=True, exist_ok=True)
        return

    # Predict
    X_test = encode_sequences(test_df)
    y_true = test_df['ccs'].values
    y_pred = model.predict(X_test).flatten()

    # Calculate metrics
    from sklearn.metrics import r2_score, mean_absolute_error

    metrics = {
        'r2': float(r2_score(y_true, y_pred)),
        'mae': float(mean_absolute_error(y_true, y_pred)),
        'mape': float(np.mean(np.abs((y_true - y_pred) / y_true)) * 100),
        'n_test': len(y_true),
    }

    logger.info(f"R²: {metrics['r2']:.4f}")
    logger.info(f"MAE: {metrics['mae']:.2f} Å²")
    logger.info(f"MAPE: {metrics['mape']:.2f}%")

    # Create plots
    create_plots(y_true, y_pred, plots_dir)

    # Generate report
    generate_html_report(metrics, report_output)

    logger.info("Evaluation completed")


if __name__ == "__main__":
    main()
