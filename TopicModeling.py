import pandas as pd
import numpy as np
from umap import UMAP
from hdbscan import HDBSCAN
from sklearn.cluster import KMeans
from sklearn.metrics import mean_squared_error
import optuna

class TextClustering:
    def __init__(self, mapper_model: UMAP, cluster_model: HDBSCAN | KMeans) -> None:
        """
        Performs clustering of texts on reduced embeddings dimension
        Args:
            mapper_model: UMAP, PCA or tSNE model for dimensionality reduction
            cluster_model: HDBSCAN or kMeans clustering model
        """
        self.mapper_model = mapper_model
        self.cluster_model = cluster_model
        self.reduced_embeddings : np.ndarray | None = None

    def set_mapper(self, mapper_model: UMAP) -> None:
        self.mapper_model = mapper_model

    def _reduce_embeddings(self, text_embeddings: np.ndarray):
        print(f"Reducing dimensionality of embeddings with {self.mapper_model} model")
        return self.mapper_model.fit_transform(text_embeddings)

    def cluster_texts(self, text_embeddings: np.ndarray):
        # Run dimensionality reduction first
        self.reduced_embeddings = self._reduce_embeddings(text_embeddings)
        print(f"Reduced Embeddings shape: {self.reduced_embeddings.shape}")

        # Cluster the reduced embeddings
        print(f"Running clustering algorithm on the data using {self.cluster_model} model")
        self.cluster_model = self.cluster_model.fit(self.reduced_embeddings)

        # Labels of each point
        clusters = self.cluster_model.labels_
        print(f"Number of clusters: {len(set(clusters))}")

        return clusters
    
    def optimize_mapper(self, text_embeddings: np.ndarray, n_trials: int = 20):
        """
        Optimizes UMAP parameters using Optuna by minimizing reconstruction error.
        Args:
            text_embeddings: High-dimensional text embeddings.
            n_trials: Number of trials for Optuna optimization.
        Returns:
            Best parameters found for UMAP.
        """
        def objective(trial):
            # Suggest hyperparameters for UMAP
            n_neighbors = trial.suggest_int("n_neighbors", 5, 50)
            min_dist = trial.suggest_float("min_dist", 0.0, 0.99)
            n_components = trial.suggest_int("n_components", 2, 50)

            # Create a new UMAP model with suggested parameters
            umap_model = UMAP(n_neighbors=n_neighbors, min_dist=min_dist, n_components=n_components, random_state=42)

            # Reduce embeddings
            reduced_embeddings = umap_model.fit_transform(text_embeddings)

            # Attempt to reconstruct the original embeddings
            try:
                reconstructed_embeddings = umap_model.inverse_transform(reduced_embeddings)
                mse = mean_squared_error(text_embeddings, reconstructed_embeddings)
            except ValueError:
                # If inverse transform fails (e.g., due to small n_components), return a high MSE
                mse = float("inf")

            return mse

        # Use Optuna to minimize reconstruction error
        sampler = optuna.samplers.TPESampler(n_startup_trials=5, seed=42)
        study = optuna.create_study(direction="minimize", sampler=sampler)
        study.optimize(objective, n_trials=n_trials)

        # Best parameters
        best_params = study.best_params
        print(f"Best UMAP parameters: {best_params}")

        # Update the mapper model with best parameters
        self.mapper_model = UMAP(**best_params, random_state=42)

        return best_params
