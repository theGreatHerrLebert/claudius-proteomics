"""
Rules for model training via imspy-predictors.

Models are trained on database snapshots for reproducibility.
"""


rule train_ccs_model:
    """
    Train CCS prediction model on a database snapshot.

    Uses the ionmob architecture (GRU-based) for CCS prediction.
    """
    input:
        data=lambda wildcards: f"database/snapshots/{config.get('snapshot', 'v1.0')}/training_set.parquet"
    output:
        model="models/{model_name}/model.h5",
        metrics="models/{model_name}/metrics.json",
        config="models/{model_name}/config.yaml"
    params:
        test_split=config["training"]["test_split"],
        val_split=config["training"]["validation_split"],
        batch_size=config["training"]["batch_size"],
        epochs=config["training"]["epochs"],
        patience=config["training"]["early_stopping_patience"],
        snapshot=lambda wildcards: config.get('snapshot', 'v1.0')
    singularity:
        config["containers"]["imspy"]
    resources:
        mem_mb=32000,
        time="4:00:00",
        gpu=1,
        partition="gpu"
    log:
        "logs/train/{model_name}.log"
    script:
        "../scripts/train_ccs.py"


rule evaluate_trained_model:
    """
    Evaluate trained model on held-out test set.

    Produces detailed metrics and visualization plots.
    """
    input:
        model="models/{model_name}/model.h5",
        data=lambda wildcards: f"database/snapshots/{config.get('snapshot', 'v1.0')}/training_set.parquet"
    output:
        report="models/{model_name}/evaluation_report.html",
        plots=directory("models/{model_name}/plots")
    params:
        snapshot=lambda wildcards: config.get('snapshot', 'v1.0')
    singularity:
        config["containers"]["imspy"]
    resources:
        mem_mb=16000,
        time="1:00:00",
        gpu=1
    log:
        "logs/evaluate/{model_name}.log"
    script:
        "../scripts/evaluate_model.py"


rule export_model:
    """
    Export trained model for deployment/sharing.

    Creates a standalone package with model, config, and metadata.
    """
    input:
        model="models/{model_name}/model.h5",
        metrics="models/{model_name}/metrics.json",
        config="models/{model_name}/config.yaml"
    output:
        package="models/{model_name}/export/{model_name}_export.tar.gz"
    run:
        import tarfile
        from pathlib import Path

        export_dir = Path(output.package).parent
        export_dir.mkdir(parents=True, exist_ok=True)

        with tarfile.open(output.package, "w:gz") as tar:
            tar.add(input.model, arcname=f"{wildcards.model_name}/model.h5")
            tar.add(input.metrics, arcname=f"{wildcards.model_name}/metrics.json")
            tar.add(input.config, arcname=f"{wildcards.model_name}/config.yaml")

        print(f"Model exported to {output.package}")
