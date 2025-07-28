import logging
from pathlib import Path
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.metrics.pairwise import cosine_similarity
from scipy.stats import pearsonr
from torch_geometric.loader import DataLoader as PyGDataLoader
from mltutils.models import LatentGeneExpressionGNN, GeneExpressionAutoencoder
from mltutils.utils import (
    MoleculeGeneLatentDataset,
    load_gene_expression_data,
    train_autoencoder,
    train,
    predict_gene_expression,
)

logging.getLogger("deepchem").setLevel(logging.WARNING)

cell_line_id = "JURKAT"

# Configuration Section
CONFIG = {
    "data_folder": Path("data-mlt/"),
    "train_autoencoder": True,  # Whether to train the autoencoder from scratch
    "train_prediction_model": True,  # Whether to train the prediction model from scratch
    "autoencoder_weights_path": f"autoencoder_weights_{cell_line_id}.pth",
    "prediction_model_weights_path": f"prediction_model_weights_{cell_line_id}.pth",
    "learning_rate": 3e-4,
    "weight_decay": 1e-6,
    "num_epochs": 30,
    "hidden_dim": 128,
    "cell_line_embedding_dim": 8,
    "batch_size": 32,
    "batch_size_ae": 32,
    "latent_dim": 64,
    "device": torch.device("cuda" if torch.cuda.is_available() else "cpu"),
    "n_epochs_autoencoder": 30,
    "autoencoder_lr": 1e-3,
    "autoencoder_weight_decay": 1e-6,
    "train_smiles_file": "MLT.csv",
    "test_smiles_file": "test.csv",
    "test_smiles_cell_id": cell_line_id,
    "disease_data_file": "disease.csv",
    "output_file": None,
}


# Setup Logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(message)s")
logger = logging.getLogger(__name__)


# Function Definitions
def load_data(config):
    """Loads data and splits into training and test sets."""
    logger.info("Loading data...")
    data_path = config["data_folder"] / config["train_smiles_file"]
    return load_gene_expression_data(data_path)


def preprocess_data(expr_train, expr_test):
    """Normalizes data using StandardScaler."""
    logger.info("Normalizing gene expression data...")
    scaler = StandardScaler()
    expr_train_normalized = scaler.fit_transform(expr_train)
    expr_test_normalized = scaler.transform(expr_test)

    lower_bound = np.percentile(expr_train_normalized, 1)
    upper_bound = np.percentile(expr_train_normalized, 99)
    expr_train_capped = np.clip(expr_train_normalized, lower_bound, upper_bound)
    expr_test_capped = np.clip(expr_test_normalized, lower_bound, upper_bound)

    return expr_train_capped, expr_test_capped, scaler


def train_autoencoder_model(expr_train_normalized, config):
    """Trains the autoencoder model."""
    expr_train_ae, expr_val_ae = train_test_split(
        expr_train_normalized, test_size=0.2, random_state=42
    )

    expr_train_tensor = torch.FloatTensor(expr_train_ae)
    expr_val_tensor = torch.FloatTensor(expr_val_ae)
    train_loader_ae = DataLoader(
        expr_train_tensor, batch_size=config["batch_size_ae"], shuffle=True
    )
    val_loader_ae = DataLoader(
        expr_val_tensor, batch_size=config["batch_size_ae"], shuffle=False
    )

    autoencoder = GeneExpressionAutoencoder(
        input_dim=expr_train_normalized.shape[1], latent_dim=config["latent_dim"]
    ).to(config["device"])

    logger.info("Training autoencoder...")
    criterion = nn.MSELoss()
    optimizer = optim.Adam(
        autoencoder.parameters(),
        lr=config["autoencoder_lr"],
        weight_decay=config["autoencoder_weight_decay"],
    )
    scheduler = optim.lr_scheduler.ReduceLROnPlateau(
        optimizer, mode="min", factor=0.2, patience=5
    )

    trained_autoencoder = train_autoencoder(
        autoencoder,
        train_loader_ae,
        val_loader_ae,
        optimizer,
        criterion,
        config["device"],
        scheduler,
        config["n_epochs_autoencoder"],
    )
    return trained_autoencoder


# Autoencoder Training or Loading
def train_or_load_autoencoder(expr_train_normalized, config):
    """Trains the autoencoder or loads pre-trained weights."""
    autoencoder = GeneExpressionAutoencoder(
        input_dim=expr_train_normalized.shape[1], latent_dim=config["latent_dim"]
    ).to(config["device"])

    if config["train_autoencoder"]:
        logger.info("Training autoencoder from scratch...")
        autoencoder = train_autoencoder_model(expr_train_normalized, config)
        torch.save(autoencoder.state_dict(), config["autoencoder_weights_path"])
    else:
        logger.info(
            f"Loading autoencoder weights from {config['autoencoder_weights_path']}..."
        )
        try:
            autoencoder.load_state_dict(torch.load(config["autoencoder_weights_path"]))
            autoencoder.eval()
        except FileNotFoundError:
            logger.error(
                "Pre-trained autoencoder weights not found. Training from scratch..."
            )
            autoencoder = train_autoencoder_model(expr_train_normalized, config)
            torch.save(autoencoder.state_dict(), config["autoencoder_weights_path"])

    return autoencoder


# Prediction Model Training or Loading
def train_or_load_prediction_model(
    train_loader, test_loader, autoencoder, scaler, num_cell_lines, config
):
    """Trains the prediction model or loads pre-trained weights."""
    num_node_features = train_loader.dataset[0][0].num_node_features

    model = LatentGeneExpressionGNN(
        num_node_features=num_node_features,
        hidden_dim=config["hidden_dim"],
        num_cell_lines=num_cell_lines,
        cell_line_embedding_dim=config["cell_line_embedding_dim"],
        latent_dim=config["latent_dim"],
    ).to(config["device"])

    if config["train_prediction_model"]:
        logger.info("Training prediction model from scratch...")
        model, final_metrics = train_prediction_model(
            train_loader, test_loader, autoencoder, scaler, num_cell_lines, config
        )
        torch.save(model.state_dict(), config["prediction_model_weights_path"])
    else:
        logger.info(
            f"Loading prediction model weights from {config['prediction_model_weights_path']}..."
        )
        try:
            model.load_state_dict(torch.load(config["prediction_model_weights_path"]))
            model.eval()
            final_metrics = None
        except FileNotFoundError:
            logger.error(
                "Pre-trained prediction model weights not found. Training from scratch..."
            )
            model, final_metrics = train_prediction_model(
                train_loader, test_loader, autoencoder, scaler, num_cell_lines, config
            )
            torch.save(model.state_dict(), config["prediction_model_weights_path"])

    return model, final_metrics


def encode_data(autoencoder, data_normalized, config):
    """Encodes data using the trained autoencoder."""
    logger.info("Encoding data using the trained autoencoder...")
    autoencoder.eval()
    with torch.no_grad():
        encoded_data = (
            autoencoder.encoder(torch.FloatTensor(data_normalized).to(config["device"]))
            .detach()
            .cpu()
            .numpy()
        )
    return encoded_data


def train_prediction_model(
    train_loader, test_loader, autoencoder, scaler, num_cell_lines, config
):
    """Trains the graph-based gene expression prediction model."""
    logger.info("Training prediction model...")
    num_node_features = train_loader.dataset[0][0].num_node_features

    model = LatentGeneExpressionGNN(
        num_node_features=num_node_features,
        hidden_dim=config["hidden_dim"],
        num_cell_lines=num_cell_lines,
        cell_line_embedding_dim=config["cell_line_embedding_dim"],
        latent_dim=config["latent_dim"],
    ).to(config["device"])

    criterion = nn.MSELoss()
    optimizer = optim.Adam(
        model.parameters(),
        lr=config["learning_rate"],
        weight_decay=config["weight_decay"],
    )
    scheduler = optim.lr_scheduler.CosineAnnealingLR(
        optimizer, T_max=config["num_epochs"]
    )

    final_metrics = train(
        model=model,
        train_loader=train_loader,
        test_loader=test_loader,
        criterion=criterion,
        optimizer=optimizer,
        scheduler=scheduler,
        num_epochs=config["num_epochs"],
        device=config["device"],
        autoencoder=autoencoder,
        scaler=scaler,
    )
    return model, final_metrics


def main(config):
    # Load and preprocess data
    (
        graph_data,
        cell_lines,
        gene_expression,
        num_cell_lines,
        gene_names,
        cell_line_encoder,
    ) = load_data(config)
    (
        graphs_train,
        graphs_test,
        cell_lines_train,
        cell_lines_test,
        expr_train,
        expr_test,
    ) = train_test_split(
        graph_data, cell_lines, gene_expression, test_size=0.2, random_state=42
    )

    expr_train_normalized, expr_test_normalized, scaler = preprocess_data(
        expr_train, expr_test
    )
    # Train autoencoder
    autoencoder = train_or_load_autoencoder(expr_train_normalized, config)

    # Encode data
    z_train = encode_data(autoencoder, expr_train_normalized, config)
    z_test = encode_data(autoencoder, expr_test_normalized, config)

    # Create datasets and data loaders
    train_dataset = MoleculeGeneLatentDataset(graphs_train, cell_lines_train, z_train)
    test_dataset = MoleculeGeneLatentDataset(graphs_test, cell_lines_test, z_test)

    train_loader = PyGDataLoader(
        train_dataset, batch_size=config["batch_size"], shuffle=True
    )
    test_loader = PyGDataLoader(test_dataset, batch_size=config["batch_size"])

    # Train prediction model
    model, final_metrics = train_or_load_prediction_model(
        train_loader, test_loader, autoencoder, scaler, num_cell_lines, config
    )

    # Test and evaluate
    logger.info("Making predictions and evaluating...")
    df_test = pd.read_csv(config["data_folder"] / config["test_smiles_file"])
    smiles_cell_line_pairs = [(smile, config["test_smiles_cell_id"]) for smile in df_test["SMILES"].values]
    predictions = predict_gene_expression(
        model=model,
        autoencoder=autoencoder,
        smiles_cell_line_pairs=smiles_cell_line_pairs,
        cell_line_encoder=cell_line_encoder,
        device=config["device"],
        scaler=scaler,
    )

    assert (
        len(gene_names) == predictions.shape[1]
    ), "Gene names must match the number of predicted genes."

    # Compare predictions with disease data
    df_disease = pd.read_csv(config["data_folder"] / config["disease_data_file"])
    disease = df_disease["exprs"].values

    # Find common genes
    common_genes = list(set(gene_names).intersection(set(df_disease["Gene"].values)))
    predictions_common = predictions[:, [gene_names.index(gene) for gene in common_genes]]

    scaler = StandardScaler()
    predictions_normalized = scaler.fit_transform(predictions_common)
    disease_normalized = scaler.transform(disease.reshape(1, -1)).flatten()

    similarity_scores = cosine_similarity(
        predictions_normalized, disease_normalized.reshape(1, -1)
    ).flatten()
    opposite_cosine_scores = -similarity_scores

    opposite_pearson_scores = np.array(
        [
            -pearsonr(predictions_normalized[i], disease_normalized)[0]
            for i in range(predictions_normalized.shape[0])
        ]
    )
    opposite_pearson_scores = np.nan_to_num(opposite_pearson_scores)

    # Create and sort results
    df_result = pd.DataFrame(
        {
            "cid": df_test["cid"],
            "Name": df_test["cmpdname"],
            "SMILES": df_test["SMILES"],
            "Negative_Cosine_Similarity": opposite_cosine_scores,
            "Negative_Pearson_Correlation": opposite_pearson_scores,
            "Average_Score": opposite_cosine_scores + opposite_pearson_scores,
        }
    )
    df_result.sort_values("Average_Score", ascending=True, inplace=True)
    df_result.reset_index(drop=True, inplace=True)

    gene_expression_df = pd.DataFrame(data=predictions_common, columns=common_genes)
    df_result = pd.concat([df_result, gene_expression_df], axis=1)
    return df_result, autoencoder, model, final_metrics


if __name__ == "__main__":
    df_result, autoencoder, model, final_metrics = main(CONFIG)
    # print(df_result)
    df_result.to_csv(CONFIG["data_folder"] / f"MLT_result_{CONFIG['test_smiles_cell_id']}.csv", index=False)
