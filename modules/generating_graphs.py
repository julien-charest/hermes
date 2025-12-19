#!/usr/bin/env python3

####################################################################
# Hermes v1.2 - Open-source mining tool for open-access literature #
# 2025-12-18                                                       #
# Written by Julien Charest & Katarina Priselac                    #
####################################################################

import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.ticker import MaxNLocator
from itertools import cycle
import pandas as pd
import numpy as np


def graph_pub_years(dataframe, request_id, dir):

    """Generate and save a line plot showing how many retrieved articles were published per year."""

    # Define color palette for the figure
    palette = "#0C1927"

    # Count number of articles per publication year
    years_count = dataframe.groupby(["Year"])["PMID"].count()

    # Create the figure
    plt.figure(figsize = (10, 4))
    plt.plot(years_count, color = palette)

    # Configure visual elements
    plt.title("Articles Retrieved per Publication Year")
    plt.xlabel("Publication Year")
    ax = plt.gca()
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    plt.ylabel("Number of Articles")

    # Improve layout and save the final output
    plt.tight_layout()
    plt.savefig("{dir}/{request}_pubyears.png".format(dir = dir, request = request_id))
    plt.close()


def graph_associated_keywords(dataframe, request_id, dir):

    """Generate and save a bar chart showing the top 25 most frequently occurring associated keywords found in the dataset."""

    # Define color palette for the figure
    palette = ["#0C1927", "#152C45", "#1F4063", "#4B6682", "#788CA1", "#A5B2C0"]
    colors = cycle(palette)

    # Extract keyword lists from dataframe
    keywords = dataframe["Associated Keywords"]
    
    # Flatten list of keyword lists into a single sequence
    keywords_list = []
    for row in keywords:
        for i in row:
            keywords_list.append(i)
    
    # Create a dataframe of keywords for aggregation
    keyword_df = pd.DataFrame(np.array([[i, "NA"] for i in keywords_list]), columns= ["Keywords", "Count"])

    # Count frequency of each unique keyword and sort by occurrence
    keyword_df = keyword_df.groupby(["Keywords"]).count().sort_values(by= "Count", ascending= False)

    # Select top 25 most frequently occurring keywords
    top_25 = keyword_df.head(25)

    # Create the figure
    plt.figure(figsize = (10, 6))
    plt.bar(top_25.index.tolist(), top_25["Count"].squeeze(), color = [next(colors) for i in range(25)])
    
    # Configure visual elements
    plt.title("Top 25 Keywords Associated with Retrieved Articles")
    ax = plt.gca()
    ax.tick_params(axis = "x", labelsize = 8)
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    plt.xticks(rotation = -90)
    plt.ylabel("Number of Articles")
    plt.xlabel("Keywords")

    # Improve layout and save the final output
    plt.tight_layout()
    plt.savefig("{dir}/{request}_asskeywords.png".format(dir = dir, request = request_id))
    plt.close()

def graph_stats_summary(dataframe, request_id, dir):

    """
    Generate and save a summary figure showing statistics for:
        - Citations per article
        - Keyword counts per keyword per article
        - Scores per article

    The figure contains three subplots with boxplots and overlaid scatter points.

    """

    # Define color palette for the figure
    palette = ["#0C1927", "#152C45", "#1F4063", "#4B6682", "#788CA1", "#A5B2C0"]
    colors = cycle(palette)

    # Create overall figure and grid layout
    plt.figure(figsize = (10, 5))
    plt.suptitle("Hermes Report Summary")
    gs = gridspec.GridSpec(1, 3, width_ratios = [1, 4, 1], wspace=0.5)

    # --- Subplot 1: citations statistics ---
    fig1 = plt.subplot(gs[0])
    labels = ["Citations"]

    # Prepare data for boxplot and scatter
    y = [dataframe["Citations"].squeeze().to_list()]
    x = []
    for i in range(len(y)):
        x.append([i + 1 for z in range(len(y[0]))])

    # Create the subplot
    fig1.boxplot(y, labels = labels)
    fig1.scatter(x, y, color = "#0C1927", alpha=0.5)

    # Configure visual elements
    fig1.set_xticks([])
    fig1.tick_params(axis = "y", labelsize = 8)
    fig1.set_xlabel("Citations")
    fig1.set_ylabel("Number of Citations")

    # --- Subplot 2: keyword count statistics ---
    fig2 = plt.subplot(gs[1])

    # Extract per-article keyword count dictionaries
    keywords_count_list = dataframe["Keywords Count"].squeeze().to_list()

    # Use keys from the first row as keyword labels
    labels = [i for i in keywords_count_list[0].keys()]
    
    # Build list of values per keyword across all articles
    y = []
    for i in labels:
        y.append([z[i] for z in keywords_count_list])

    # Flatten x and y for scatter plotting
    x = []
    for i in range(len(y)):
        x.extend([i + 1] * len(y[i]))
    y_flat = [val for sublist in y for val in sublist]

    # Assign colors per keyword group
    c = []
    for i in range(len(y)):
        c.extend([next(colors)] * len(y[i]))

    # Create the subplot
    fig2.boxplot(y, labels = labels)
    fig2.scatter(x, y_flat, color = c, alpha=0.5)

    # Configure visual elements
    fig2.tick_params(axis = "x", labelsize = 8, labelrotation = 0)
    fig2.tick_params(axis = "y", labelsize = 8)
    fig2.set_xlabel("Keywords")
    fig2.set_ylabel("Keyword Count per Article")

    # --- Subplot 3: score statistics ---
    fig3 = plt.subplot(gs[2])
    labels = ["Score"]

    # Prepare data for boxplot and scatter
    y = [dataframe["Score"].squeeze().to_list()]
    x = []
    for i in range(len(y)):
        x.append([i + 1 for z in range(len(y[0]))])

    # Create the subplot
    fig3.boxplot(y, labels = labels)
    fig3.scatter(x, y, color = next(colors), alpha=0.5)

    # Configure visual elements
    fig3.set_xticks([])
    fig3.tick_params(axis = "y", labelsize = 8)
    fig3.set_xlabel("Score")
    fig3.set_ylabel("Score")

    # Save final summary figure
    plt.savefig("{dir}/{request}_stats_summary.png".format(dir = dir, request = request_id))
    plt.close()