import os, sys
import pandas as pd
from sentence_transformers import SentenceTransformer
from umap import UMAP
from hdbscan import HDBSCAN
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from itertools import cycle


from DataLoader import MathProblemsDataset

def main():
    # Get data
    ds = MathProblemsDataset(dataset_name="AI-MO/NuminaMath-TIR", partition_name='train')

    dataset = ds.get_data()
    problems = dataset['problem'][:4000]

    #for ipbm, problem in enumerate(problems):
    #    print(f"Problem {ipbm}: {problem}")

    encoder_model = SentenceTransformer("multi-qa-MiniLM-L6-cos-v1") # config_kwargs={"use_memory_efficient_attention": False}
    problems_embedded = encoder_model.encode(problems, show_progress_bar=True)
    print(f"Embedding shape: {problems_embedded.shape}")

    umap_model = UMAP(n_components=2, min_dist=0.0, metric='cosine', random_state=42)
    problems_reduced_embeddings = umap_model.fit_transform(problems_embedded)
    print(f"Reduced Embeddings shape: {problems_reduced_embeddings.shape}")

    hdbscan_model = HDBSCAN(min_cluster_size=40, metric="euclidean", cluster_selection_method="eom").fit(problems_reduced_embeddings)
    clusters = hdbscan_model.labels_
    print(f"Number of clusters: {len(set(clusters))}")


    df = pd.DataFrame(problems_reduced_embeddings, columns=['x', 'y'])
    df['cluster'] = [str(cl) for cl in clusters]
    df_clusters = df.loc[df.cluster!='-1']
    df_outliers = df.loc[df.cluster=='-1']

    plt.figure()
    plt.scatter(df_outliers.x, df_outliers.y, alpha=0.5, s=20, c='grey', label='outliers')
    colors = cycle(cm.tab10.colors)
    for cl in df_clusters.cluster.sort_values().unique():
        mask = df_clusters.cluster==cl
        color = next(colors)
        plt.scatter(
            df_clusters[mask].x, df_clusters[mask].y, c=color, alpha=0.6, s=20, label='cluster '+cl
        )
    plt.legend()
    plt.savefig("1.png")


if __name__ == "__main__":
    print('Running ...')
    main()