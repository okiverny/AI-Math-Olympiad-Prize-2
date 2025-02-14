import os, sys
import pandas as pd
import numpy as np
from sentence_transformers import SentenceTransformer
from umap import UMAP
from hdbscan import HDBSCAN
from sklearn.cluster import KMeans

from DataLoader import MathProblemsDataset
#from DataEmbedding import MathDataEmbedding
from TopicModeling import TextClustering
from plotting import (
    plot_clusters,
    plot_clustering_tree,
    BERTopic_visualize_barchart
)

np.random.seed(42)

def main():
    # Get data
    mpds = MathProblemsDataset(dataset_name="AI-MO/NuminaMath-TIR", partition_name='train[:40000]')

    dataset = mpds.data
    problems = dataset['problem']

    # problems_embedded = MathDataEmbedding("multi-qa-MiniLM-L6-cos-v1").encode(problems, show_progress_bar=True)
    encoder_model = SentenceTransformer("multi-qa-MiniLM-L6-cos-v1")
    problems_embedded = encoder_model.encode(problems, show_progress_bar=True)
    print(f"Embedding shape: {problems_embedded.shape}")

    # Initialize TextClustering with a default UMAP and clustering model
    clustering_model = TextClustering(
        mapper_model = UMAP(n_components=5, min_dist=0.0, metric='cosine', random_state=42),
        cluster_model = HDBSCAN(min_cluster_size=400, metric="euclidean", cluster_selection_method="eom", prediction_data=True)
        #cluster_model = KMeans(n_clusters=4)
    )

    # Uncomment to optimize UMAP parameters (be warned this is quite expensive computationally!!)
    # best_umap_params = clustering_model.optimize_mapper(problems_embedded, n_trials=5)

    # Cluster texts using optimized UMAP
    clusters = clustering_model.cluster_texts(problems_embedded)

    mpds.data = dataset.add_column(name="clusters", column=clusters)
    #print(mpds.data.features)

    clusters, probs, keywords = clustering_model.BERTopic_train(encoder_model, problems, problems_embedded)
    mpds.data = dataset.add_column(name="clusters_BERTopic", column=clusters)

    # Exemine clusters (printing problems with or without solutions)
    mpds.exemine_clusters('clusters_BERTopic', show_solution=False)

    # Rerun the UMAP for a 2-dimensional plot
    problems_reduced_embeddings = UMAP(n_components=2, min_dist=0.0, metric='cosine', random_state=42).fit_transform(problems_embedded)
    print(f"Reduced Embeddings shape: {problems_reduced_embeddings.shape}")


    df = pd.DataFrame(mpds.data)
    print(df.head())

    # Some illustrative plots
    df = pd.DataFrame(problems_reduced_embeddings, columns=['x', 'y'])
    df['cluster'] = [str(cl) for cl in clusters]
    plot_clusters(df)
    print('Here...')
    plot_clustering_tree(clustering_model.cluster_model)
    print('Here2...')
    BERTopic_visualize_barchart(clustering_model.topic_model, top_n_topics=3)

    # My queries
    print('Playing with queries...')
    queries = ["This problem is related to geometry and involves triangles and circles.",
               "This problem is related to the analysis of functions.",
               "This problem is related to operations with complex numbers.",
               "This problem is related to series expansions",
               "This problem is related to finding numbers with certain properties of their digits.",
               "This problem is about solving polynomial equations.",
               "This problem is related to sequences such as arithmetic and geometric progressions.",
               "This problem is related to sets and subsets.",
               "This problem is related to inequalities.",
               "This problem involves recurrence relations.",
               "This problem is related to combinatorial algebra.",
               "This problem involves trigonometric functions.",
               "This problem is related to Graph Theory.",
               "This problem is related to Probability and Expected Value.",
               "This problem is involves real-world situations and requires analytical thinking."]


    query_topics, query_probs = clustering_model.topic_model.transform(queries)
    for i, query in enumerate(queries):
        print(f"\nQuery: {query}")
        print(f"Class: {query_topics[i]}, Prob: {query_probs[i]}")
        print(f"Keywords: {clustering_model.topic_model.get_topic(query_topics[i])}")




if __name__ == "__main__":
    print('Running ...')
    main()