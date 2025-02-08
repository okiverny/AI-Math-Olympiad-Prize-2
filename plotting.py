import matplotlib.pyplot as plt
import matplotlib.cm as cm
from itertools import cycle
import pandas as pd
from bertopic import BERTopic

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=UserWarning)

def plot_clusters(df: pd.DataFrame) -> None:
    """
    Function to plot clusters on 2-dimensional projection from clustering algorithm.
    Args:
        df (pd.DataFrame): dataframe with `x` and `y` columns as well as with `cluster` columns representing cluster labels
    """
    df_clusters = df.loc[df.cluster!='-1']
    df_outliers = df.loc[df.cluster=='-1']
    print(f"Number of outliers: {len(df_outliers)}")

    plt.figure()
    if len(df_outliers)>0:
        plt.scatter(df_outliers.x, df_outliers.y, alpha=0.5, s=20, c='grey', label='outliers')

    colors = cycle(cm.tab10.colors) if len(df.cluster.unique())<10 else cycle(cm.tab20.colors)
    for cl in df_clusters.cluster.sort_values().unique():
        mask = df_clusters.cluster==cl
        color = next(colors)
        plt.scatter(
            df_clusters[mask].x, df_clusters[mask].y, color=color, alpha=0.6, s=20, label='cluster '+cl
        )
    plt.legend()
    plt.savefig("clusters.png")

def plot_clustering_tree(clustering_model) -> None:
    plt.figure()
    clustering_model.condensed_tree_.plot(select_clusters=True)
    plt.savefig("clusters_tree.png")

def BERTopic_visualize_barchart(
        topic_model: BERTopic,
        top_n_topics: int = 3,
        n_words: int = 7
) -> None:
    fig = topic_model.visualize_barchart(top_n_topics=top_n_topics, n_words=n_words, height=700)
    fig.write_image("topics.png")