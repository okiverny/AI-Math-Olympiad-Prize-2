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

    # Initialize TextClustering with a default UMAP and clustering model
    clustering_model = TextClustering(
        mapper_model = UMAP(n_components=5, min_dist=0.0, metric='cosine', random_state=42),
        #cluster_model = HDBSCAN(min_cluster_size=40, metric="euclidean", cluster_selection_method="eom")
        cluster_model = KMeans(n_clusters=4)
    )

    # Uncomment to optimize UMAP parameters (be warned this is quite expensive computationally!!)
    # best_umap_params = clustering_model.optimize_mapper(problems_embedded, n_trials=5)

    # Cluster texts using optimized UMAP
    clusters = clustering_model.cluster_texts(problems_embedded)

    mpds.data = dataset.add_column(name="clusters", column=clusters)
    print(mpds.data.features)

    # Exemine clusters
    for cluster in set(clusters):  # list(set(clusters)).sort()
        print(f"======= Cluster {cluster} ==============")
        for problem in mpds.data.filter(lambda example: example["clusters"]==cluster)['problem'][:5]:
            print(problem)


    # Rerun the UMAP for a 2-dimensional plot
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