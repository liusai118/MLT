import torch
import torch.nn as nn
from torch_geometric.nn import GCNConv, global_add_pool


class GeneExpressionAutoencoder(nn.Module):
    def __init__(self, input_dim=978, latent_dim=50):
        super().__init__()

        self.encoder = nn.Sequential(
            nn.Linear(input_dim, 256),
            nn.ReLU(),
            nn.Dropout(0.2),
            nn.Linear(256, latent_dim),
        )

        self.decoder = nn.Sequential(
            nn.Linear(latent_dim, 256),
            nn.ReLU(),
            nn.Dropout(0.2),
            nn.Linear(256, input_dim),
        )

        self._initialize_weights()

    def _initialize_weights(self):
        for module in self.modules():
            if isinstance(module, nn.Linear):
                nn.init.kaiming_normal_(module.weight, nonlinearity="relu")
                if module.bias is not None:
                    nn.init.zeros_(module.bias)

    def forward(self, x):
        z = self.encoder(x)
        reconstructed_x = self.decoder(z)
        return reconstructed_x


class LatentGeneExpressionGNN(nn.Module):
    def __init__(
        self,
        num_node_features,
        hidden_dim,
        num_cell_lines,
        cell_line_embedding_dim,
        latent_dim,
    ):
        super().__init__()

        self.gnn1 = GCNConv(num_node_features, hidden_dim)
        self.gnn2 = GCNConv(hidden_dim, hidden_dim)

        self.cell_line_embedding = nn.Embedding(num_cell_lines, cell_line_embedding_dim)

        self.combine_layer = nn.Linear(hidden_dim + cell_line_embedding_dim, hidden_dim)
        self.ln_combine = nn.LayerNorm(hidden_dim)

        self.fc1 = nn.Linear(hidden_dim, hidden_dim)
        self.ln_fc1 = nn.LayerNorm(hidden_dim)
        self.fc2 = nn.Linear(hidden_dim, latent_dim)

        self.dropout = nn.Dropout(0.2)
        self.relu = nn.ReLU()

        self.reset_parameters()

    def reset_parameters(self):
        for module in self.modules():
            if isinstance(module, nn.Linear):
                nn.init.kaiming_uniform_(module.weight)
                if module.bias is not None:
                    nn.init.zeros_(module.bias)
            elif isinstance(module, nn.Embedding):
                nn.init.normal_(module.weight, mean=0, std=1)
            elif isinstance(module, nn.LayerNorm):
                nn.init.ones_(module.weight)
                nn.init.zeros_(module.bias)

    def forward(self, data, cell_lines):
        x, edge_index, batch = data.x, data.edge_index, data.batch
        x = self.gnn1(x, edge_index).relu()
        x = self.gnn2(x, edge_index).relu()

        graph_embedding = global_add_pool(x, batch)

        cell_line_embedded = self.cell_line_embedding(cell_lines)

        combined_input = torch.cat([graph_embedding, cell_line_embedded], dim=1)
        combined_encoded = self.relu(
            self.ln_combine(self.combine_layer(combined_input))
        )

        out = self.dropout(self.relu(self.fc1(combined_encoded)))
        out = self.ln_fc1(out)
        latent_pred = self.fc2(out)
        return latent_pred
