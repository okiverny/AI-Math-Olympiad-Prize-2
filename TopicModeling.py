import pandas as pd
import numpy as np
from umap import UMAP
from hdbscan import HDBSCAN
from sklearn.cluster import KMeans
from bertopic import BERTopic
from sentence_transformers import SentenceTransformer
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.metrics import mean_squared_error
import optuna
from typing import List, Tuple
from math_vocab import math_words

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
        self.topic_model: BERTopic | None = None
        self.reduced_embeddings: np.ndarray | None = None

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
            n_neighbors = trial.suggest_int("n_neighbors", 40, 400)
            min_dist = trial.suggest_float("min_dist", 0.0, 0.1)
            n_components = trial.suggest_int("n_components", 5, 10)

            # Create a new UMAP model with suggested parameters
            umap_model = UMAP(n_neighbors=n_neighbors, min_dist=min_dist, n_components=n_components, metric='cosine', random_state=42)

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
        trial = study.best_trial
        best_params = study.best_params
        print(f"Best UMAP parameters: {best_params}")

        print('Number of finished trials: ', len(study.trials))
        print('Best trial:')
        print('   Trial id:', trial.number)
        print('   Score:', trial.value)
        print('Params:')

        for key, value in trial.params.items():
            print('   {}: {}'.format(key, value))

        # Update the mapper model with best parameters
        self.mapper_model = UMAP(**best_params, metric='cosine', random_state=42)

        return best_params

    def BERTopic_train(self, embedding_model: SentenceTransformer, documents: List[str], text_embeddings: np.ndarray) -> Tuple[List[int], np.ndarray]:

        # Create CountVercotizer  and remove stop words
        #vectorizer_model = CountVectorizer(ngram_range=(1, 2), stop_words="english")
        vectorizer_model = CountVectorizer(ngram_range=(1, 4), vocabulary=math_words)

        # Train topic model with our previously defined mapper and clustering models models
        self.topic_model = BERTopic(
            embedding_model=embedding_model,
            umap_model=self.mapper_model,
            hdbscan_model=self.cluster_model,
            vectorizer_model=vectorizer_model,
            verbose=True
        )
        
        topics, probs = self.topic_model.fit_transform(documents, text_embeddings)
        #print(topics)
        print(probs)

        #print(self.topic_model.get_topic(0))
        #print(self.topic_model.get_topic(1))
        #print(self.topic_model.get_topic(2))

        df_topics = self.topic_model.get_topic_info()
        print(df_topics)
        print(df_topics.Representative_Docs[1])

        return topics, probs
