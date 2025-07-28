import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import torch
from torch.utils.data import Dataset
from torch_geometric.data import Data
from torch_geometric.loader import DataLoader as PyGDataLoader
import deepchem as dc
from scipy.stats import pearsonr
from sklearn.metrics import r2_score


class MoleculeGeneLatentDataset(Dataset):
    def __init__(self, graphs, cell_lines, latent_vectors):
        self.graphs = graphs
        self.cell_lines = torch.LongTensor(cell_lines)
        self.latent_vectors = torch.FloatTensor(latent_vectors)

    def __len__(self):
        return len(self.graphs)

    def __getitem__(self, idx):
        return self.graphs[idx], self.cell_lines[idx], self.latent_vectors[idx]


def mol_to_pyg_graph(mol):
    node_features = torch.tensor(mol.node_features, dtype=torch.float)
    edge_index = torch.tensor(mol.edge_index, dtype=torch.long)
    return Data(x=node_features, edge_index=edge_index)


def load_gene_expression_data(smiles_expr_file, input_cols=None):
    data = pd.read_csv(smiles_expr_file)
    smiles = data["SMILES"].values
    cell_lines = data["cell_id"].values
    if not input_cols:
        input_cols = ["SMILES", "cell_id"]
    gene_expression = data.drop(input_cols, axis=1).values
    gene_names = data.columns[len(input_cols):].tolist()

    featurizer = dc.feat.MolGraphConvFeaturizer()
    graphs = featurizer.featurize(smiles)

    graph_data = [mol_to_pyg_graph(mol) for mol in graphs]

    cell_line_encoder = {line: idx for idx, line in enumerate(np.unique(cell_lines))}
    encoded_cell_lines = np.array([cell_line_encoder[line] for line in cell_lines])

    return (
        graph_data,
        encoded_cell_lines,
        gene_expression,
        len(cell_line_encoder),
        gene_names,
        cell_line_encoder,
    )


def train_one_epoch(
    model, dataloader, criterion, optimizer, device, autoencoder, scaler
):
    model.train()
    losses = []
    all_preds = []
    all_targets = []

    for batch in dataloader:
        data, cell_lines, latent_vectors = batch
        data = data.to(device)
        cell_lines = cell_lines.to(device)
        latent_vectors = latent_vectors.to(device)

        optimizer.zero_grad()

        latent_pred = model(data, cell_lines)
        loss = criterion(latent_pred, latent_vectors)

        loss.backward()
        optimizer.step()

        losses.append(loss.item())

        reconstructed_pred = autoencoder.decoder(latent_pred)
        reconstructed_true = autoencoder.decoder(latent_vectors)

        all_preds.append(reconstructed_pred.detach().cpu().numpy())
        all_targets.append(reconstructed_true.detach().cpu().numpy())

    avg_loss = np.mean(losses)

    all_preds = scaler.inverse_transform(np.concatenate(all_preds))
    all_targets = scaler.inverse_transform(np.concatenate(all_targets))
    r2 = r2_score(all_targets, all_preds)
    pearson_corr, _ = pearsonr(all_targets.flatten(), all_preds.flatten())

    return avg_loss, r2, pearson_corr


def evaluate(model, dataloader, criterion, device, autoencoder, scaler):
    model.eval()
    losses = []
    all_preds = []
    all_targets = []

    with torch.no_grad():
        for batch in dataloader:
            data, cell_lines, latent_vectors = batch
            data = data.to(device)
            cell_lines = cell_lines.to(device)
            latent_vectors = latent_vectors.to(device)

            latent_pred = model(data, cell_lines)
            loss = criterion(latent_pred, latent_vectors)

            losses.append(loss.item())

            reconstructed_pred = autoencoder.decoder(latent_pred)
            reconstructed_true = autoencoder.decoder(latent_vectors)

            all_preds.append(reconstructed_pred.detach().cpu().numpy())
            all_targets.append(reconstructed_true.detach().cpu().numpy())

    avg_loss = np.mean(losses)

    all_preds = scaler.inverse_transform(np.concatenate(all_preds))
    all_targets = scaler.inverse_transform(np.concatenate(all_targets))
    r2 = r2_score(all_targets, all_preds)
    pearson_corr, _ = pearsonr(all_targets.flatten(), all_preds.flatten())

    return avg_loss, r2, pearson_corr


def train(
    model,
    train_loader,
    test_loader,
    criterion,
    optimizer,
    scheduler,
    num_epochs,
    device,
    autoencoder,
    scaler,
):
    train_loss_history = []
    val_loss_history = []
    final_r2 = 0
    final_pearson = 0

    for epoch in range(num_epochs):
        train_loss, train_r2, train_pearson = train_one_epoch(
            model, train_loader, criterion, optimizer, device, autoencoder, scaler
        )
        test_loss, test_r2, test_pearson = evaluate(
            model, test_loader, criterion, device, autoencoder, scaler
        )
        scheduler.step()

        train_loss_history.append(train_loss)
        val_loss_history.append(test_loss)
        final_r2 = test_r2
        final_pearson = test_pearson

        print(
            f"Epoch {epoch + 1}/{num_epochs}, "
            f"Train Loss: {train_loss:.4f}, Train R2: {train_r2:.4f}, Train Pearson: {train_pearson:.4f}, "
            f"Test Loss: {test_loss:.4f}, Test R2: {test_r2:.4f}, Test Pearson: {test_pearson:.4f}, "
            f"Learning Rate: {scheduler.get_last_lr()[0]:.6f}"
        )
    final_metrics = {
        "Test Loss": val_loss_history[-1],
        "Test R2": final_r2,
        "Test Pearson": final_pearson,
    }

    plt.figure()
    plt.plot(range(1, num_epochs + 1), train_loss_history, label="Train Loss")
    plt.plot(range(1, num_epochs + 1), val_loss_history, label="Validation Loss")
    plt.xlabel("Epoch")
    plt.ylabel("Loss")
    plt.title("GNN Training and Validation Loss")
    plt.legend()
    plt.show()
    return final_metrics


def train_autoencoder(
    autoencoder,
    train_loader,
    val_loader,
    optimizer,
    criterion,
    device,
    scheduler,
    num_epochs,
):
    train_loss_history = []
    val_loss_history = []

    for epoch in range(num_epochs):
        autoencoder.train()
        train_losses = []
        for batch in train_loader:
            inputs = batch.to(device)
            optimizer.zero_grad()
            outputs = autoencoder(inputs)
            loss = criterion(outputs, inputs)
            loss.backward()
            optimizer.step()
            train_losses.append(loss.item())
        avg_train_loss = sum(train_losses) / len(train_losses)
        train_loss_history.append(avg_train_loss)

        autoencoder.eval()
        val_losses = []
        with torch.no_grad():
            for batch in val_loader:
                inputs = batch.to(device)
                outputs = autoencoder(inputs)
                loss = criterion(outputs, inputs)
                val_losses.append(loss.item())
        avg_val_loss = sum(val_losses) / len(val_losses)
        val_loss_history.append(avg_val_loss)

        if isinstance(scheduler, torch.optim.lr_scheduler.ReduceLROnPlateau):
            scheduler.step(avg_val_loss)
        else:
            scheduler.step()

        print(
            f"Epoch {epoch+1}/{num_epochs}, Train Loss: {avg_train_loss:.6f}, Validation Loss: {avg_val_loss:.6f}"
        )

    plt.figure()
    plt.plot(range(1, num_epochs + 1), train_loss_history, label="Train Loss")
    plt.plot(range(1, num_epochs + 1), val_loss_history, label="Validation Loss")
    plt.xlabel("Epoch")
    plt.ylabel("Loss")
    plt.title("Autoencoder Training and Validation Loss")
    plt.legend()
    plt.show()
    return autoencoder


def predict_gene_expression(
    model, autoencoder, smiles_cell_line_pairs, cell_line_encoder, device, scaler
):
    model.eval()
    autoencoder.eval()

    featurizer = dc.feat.MolGraphConvFeaturizer()

    graphs = []
    encoded_cell_lines = []
    for smile, cell_line in smiles_cell_line_pairs:
        graph = featurizer.featurize([smile])[0]
        graph = mol_to_pyg_graph(graph)
        graphs.append(graph)
        encoded_cell_lines.append(cell_line_encoder[cell_line])

    dataset = MoleculeGeneLatentDataset(
        graphs, encoded_cell_lines, np.zeros((len(graphs), model.fc2.out_features))
    )
    data_loader = PyGDataLoader(dataset, batch_size=len(graphs), shuffle=False)

    with torch.no_grad():
        for data, cell_lines, _ in data_loader:
            data = data.to(device)
            cell_lines = cell_lines.to(device)
            latent_pred = model(data, cell_lines)
            reconstructed_expr = autoencoder.decoder(latent_pred)

    gene_expr_pred = reconstructed_expr.detach().cpu().numpy()
    gene_expr_pred_original = scaler.inverse_transform(gene_expr_pred)

    return gene_expr_pred_original
