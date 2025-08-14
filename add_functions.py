#!/usr/bin/env python3

####################################################################
# Hermes v1.1 - Open-source mining tool for open-access literature #
# 2025-08-14                                                       #
# Written by Julien Charest & Katarina Priselac                    #
####################################################################

from os import path, getcwd, mkdir
from bs4 import BeautifulSoup as bs
from urllib.request import urlopen
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.ticker import MaxNLocator
from datetime import datetime
from itertools import cycle

def get_citations(pmid):
    # Returns the number of Pubmed citations (Cited by) for a Pubmed ID (pmid)
    link = "https://pubmed.ncbi.nlm.nih.gov/?linkname=pubmed_pubmed_citedin&from_uid={pmid}".format(pmid = pmid)
    try:
        with urlopen(link) as webpage:
            data = webpage.read()

        bs_data = bs(data, "lxml")
        script_data = bs_data.find_all("script")
        citations = None
        n_citations = 0
        for section in script_data:
            section_text = section.get_text()
            if "searchQuery" in section_text:
                position_query = section_text.find("searchQuery")
                position_constants = section_text.find("searchConstants")
                citations = section_text[(position_query + len('searchQuery: "')):(position_constants)].split('",', 1)[0].split(",")
                n_citations = len(citations)
        if n_citations is None:
            n_citations = 0
        return n_citations
    
    except:
        n_citations = "NaN"

def keyword_counter(keywords, article):
    keyword_count = {}
    if len(keywords) > 1:
        for keyword in keywords:
            keyword_count[keyword] = article.lower().count(keyword.lower())
    else:
        keyword = keywords[0]
        keyword_count[keyword] = article.lower().count(keyword.lower())
    return keyword_count

def create_reports_folder(working_directory, request_id, date):
    report_dir = path.join("{working_directory}".format(working_directory = working_directory), "{date}_{request_id}".format(request_id = request_id, date = date))
    figures_dir = path.join("{report_dir}".format(report_dir = report_dir), "figures")
    results_dir = path.join("{report_dir}".format(report_dir = report_dir), "results")

    if not path.exists(report_dir):
        mkdir(report_dir)
    
    if not path.exists(figures_dir):
        mkdir(figures_dir)

    if not path.exists(results_dir):
        mkdir(results_dir)

def graph_pub_years(dataframe, request_id, dir):
    years_count = dataframe.groupby(["Year"])["PMID"].count()
    plt.figure(figsize = (10, 4))
    plt.plot(years_count, color = "#0C1927")
    plt.title("Articles Retrieved per Publication Year")
    plt.xlabel("Publication Year")
    ax = plt.gca()
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    plt.ylabel("Number of Articles")
    plt.tight_layout()
    plt.savefig("{dir}/{request}_pubyears.png".format(dir = dir, request = request_id))
    plt.close()

def graph_associated_keywords(dataframe, request_id, dir):
    keywords = dataframe["Associated Keywords"]
    colors = cycle(["#0C1927", "#152C45", "#1F4063", "#4B6682", "#788CA1", "#A5B2C0"])
    keywords_list = []
    for row in keywords:
        for i in row:
            keywords_list.append(i)
    keyword_df = pd.DataFrame(np.array([[i, "NaN"] for i in keywords_list]), columns= ["Keywords", "Count"])
    keyword_df = keyword_df.groupby(["Keywords"]).count().sort_values(by= "Count", ascending= False)
    top_25 = keyword_df.head(25)
    plt.figure(figsize = (10, 6))
    plt.bar(top_25.index.tolist(), top_25["Count"].squeeze(), color = [next(colors) for i in range(25)])
    plt.title("Top 25 Keywords Associated with Retrieved Articles")
    ax = plt.gca()
    ax.tick_params(axis = "x", labelsize = 8)
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    plt.xticks(rotation = -90)
    plt.ylabel("Number of Articles")
    plt.xlabel("Keywords")
    plt.tight_layout()
    plt.savefig("{dir}/{request}_asskeywords.png".format(dir = dir, request = request_id))
    plt.close()

def graph_stats_summary(dataframe, request_id, dir):
    colors = cycle(["#0C1927", "#152C45", "#1F4063", "#4B6682", "#788CA1", "#A5B2C0"])
    plt.figure(figsize = (10, 5))
    plt.suptitle("Hermes Report Summary")
    gs = gridspec.GridSpec(1, 3, width_ratios = [1, 4, 1], wspace=0.5)

    # Create subplot for citations statistics
    fig1 = plt.subplot(gs[0])
    labels = ["Citations"]
    y = [dataframe["Citations"].squeeze().to_list()]
    x = []
    for i in range(len(y)):
        x.append([i + 1 for z in range(len(y[0]))])
    fig1.boxplot(y, labels = labels)
    fig1.scatter(x, y, color = "#0C1927", alpha=0.5)
    fig1.set_xticks([])
    fig1.tick_params(axis = "y", labelsize = 8)
    fig1.set_xlabel("Citations")
    fig1.set_ylabel("Number of Citations")

    # Create subplot for keywords count statistics
    fig2 = plt.subplot(gs[1])
    keywords_count_list = dataframe["Keywords Count"].squeeze().to_list()
    labels = [i for i in keywords_count_list[0].keys()]
    y = []
    for i in labels:
        y.append([z[i] for z in keywords_count_list])

    x = []
    for i in range(len(y)):
        x.extend([i + 1] * len(y[i]))

    y_flat = [val for sublist in y for val in sublist]

    c = []
    for i in range(len(y)):
        c.extend([next(colors)] * len(y[i]))

    fig2.boxplot(y, labels = labels)
    fig2.scatter(x, y_flat, color = c, alpha=0.5)
    fig2.tick_params(axis = "x", labelsize = 8, labelrotation = 0)
    fig2.tick_params(axis = "y", labelsize = 8)
    fig2.set_xlabel("Keywords")
    fig2.set_ylabel("Keyword Count per Article")

    # Create subplot for score statistics
    fig3 = plt.subplot(gs[2])
    labels = ["Score"]
    y = [dataframe["Score"].squeeze().to_list()]
    x = []
    for i in range(len(y)):
        x.append([i + 1 for z in range(len(y[0]))])
    fig3.boxplot(y, labels = labels)
    fig3.scatter(x, y, color = next(colors), alpha=0.5)
    fig3.set_xticks([])
    fig3.tick_params(axis = "y", labelsize = 8)
    fig3.set_xlabel("Score")
    fig3.set_ylabel("Score")
    plt.savefig("{dir}/{request}_stats_summary.png".format(dir = dir, request = request_id))
    plt.close()

def filter_string(string):
    new_string = ""
    for i in string:
        if i in "abcdefghijklmnopqrstuvwxyz" + "abcdefghijklmnopqrstuvwxyz".upper() + " 123456789.:,;.()!?=-+/%":
            new_string += i
    return str(new_string)
