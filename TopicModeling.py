import pandas as pd
import numpy as np
from umap import UMAP
from hdbscan import HDBSCAN
from typing import List


class TextClustering:
    def __init__(
        self,
        dim_red_model,
        cluster_model,
    ) -> None:
        """
        Performs clustering of texts on reduced embeddings dimension
        Args:
            dim_red_model: UMAP, PCA or tSNE model for dimensionality reduction
            cluster_model: HDBSCAN or kMeans clustering model
        """
        self.dim_red_model = dim_red_model
        self.cluster_model = cluster_model
        self.reduced_embeddings : np.ndarray

    def _reduce_embeddings(self, text_embeddings: np.ndarray):
        print(f"Reducing dimensionality of embeddings with {self.dim_red_model} model")
        return self.dim_red_model.fit_transform(text_embeddings)

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
