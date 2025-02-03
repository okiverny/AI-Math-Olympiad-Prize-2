import os, sys
import pandas as pd
from sentence_transformers import SentenceTransformer
from umap import UMAP
from hdbscan import HDBSCAN
from sklearn.cluster import KMeans

from DataLoader import MathProblemsDataset
from DataEmbedding import MathDataEmbedding
from TopicModeling import TextClustering
from plotting import plot_clusters

def main():
    # Get data
    mpds = MathProblemsDataset(dataset_name="AI-MO/NuminaMath-TIR", partition_name='train[:4000]')

    dataset = mpds.data
    problems = dataset['problem']

    # problems_embedded = MathDataEmbedding("multi-qa-MiniLM-L6-cos-v1").encode(problems, show_progress_bar=True)
    encoder_model = SentenceTransformer("multi-qa-MiniLM-L6-cos-v1")
    problems_embedded = encoder_model.encode(problems, show_progress_bar=True)
    print(f"Embedding shape: {problems_embedded.shape}")


    clustering_model = TextClustering(
        dim_red_model = UMAP(n_components=5, min_dist=0.0, metric='cosine', random_state=42),
        #cluster_model = HDBSCAN(min_cluster_size=40, metric="euclidean", cluster_selection_method="eom")
        cluster_model = KMeans(n_clusters=2)
    )
    clusters = clustering_model.cluster_texts(problems_embedded)

    mpds.data = dataset.add_column(name="clusters", column=clusters)
    print(mpds.data.features)

    # Rerun UMAP for a 2-dimensional plot
    problems_reduced_embeddings = UMAP(n_components=2, min_dist=0.0, metric='cosine', random_state=42).fit_transform(problems_embedded)
    print(f"Reduced Embeddings shape: {problems_reduced_embeddings.shape}")


    df = pd.DataFrame(mpds.data)
    print(df.head())

    df = pd.DataFrame(problems_reduced_embeddings, columns=['x', 'y'])
    df['cluster'] = [str(cl) for cl in clusters]
    plot_clusters(df)


if __name__ == "__main__":
    print('Running ...')
    main()