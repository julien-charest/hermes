#!/usr/bin/env python3

####################################################################
# Hermes v1.2 - Open-source mining tool for open-access literature #
# 2025-12-18                                                       #
# Written by Julien Charest & Katarina Priselac                    #
####################################################################

import math
from collections import Counter
from datetime import datetime

def calculate_score(row, idf_dict):
    
    """
    Calculate a weighted relevance score for an article row using:
        - TF-IDF contribution for each scoring keyword
        - Citation count (boosts score)
        - Publication age (penalizes older papers)

    Scoring function:
    score = (w_tf_idf * tfidf_sum) * (1 + w_cit * cit_count) / (1 + w_year * age)

    """
    
    # Weighting factors controlling influence of each scoring component
    w_tf_idf = 1000   # Weight on semantic relevance (TF-IDF)
    w_cit = 0.1      # Weight on boost for citation count
    w_year = 0.1     # Weight on penalty factor for older publications

    # Compute the total TF-IDF contribution for this row
    tf_dict = row['TF']
    tfidf_sum = sum(tf * idf_dict.get(term, 0.0) for term, tf in tf_dict.items())

    # Calculate publication age
    current_year = datetime.now().year
    age = current_year - row['Year']

    # Compute final score
    score = (w_tf_idf * tfidf_sum) * (1 + w_cit * row['Citations']) / (1 + w_year * age)

    return score


def scoring_articles(dataframe, keywords_mandatory):

    """ Compute TF, IDF, and final relevance scores for all articles in the dataset."""

    # Compute global IDF values across all documents
    idf = calculate_idf(dataframe)

    # Compute TF values per row
    dataframe['TF'] = dataframe.apply(calculate_tf, axis=1)

    # Compute final relevance score for each article
    dataframe['Score'] = dataframe.apply(lambda row: calculate_score(row, idf), axis=1)

    # Optional 'Mandatory Keyword' rule enforcement: documents missing any keyword receive score = 0
    if keywords_mandatory:
        dataframe['Score'] = dataframe.apply(lambda row: 0 if any(count == 0 for count in row['Keywords Count'].values()) else row['Score'], axis=1)

    return dataframe, idf


def calculate_idf(dataframe):
    
    """Compute the Inverse Document Frequency (IDF) values for scoring keywords in a dataset."""

    # Total number of documents in the dataset
    N = len(dataframe)
    
    # Counter to track how many documents contain each term
    doc_freq = Counter()
    for counts in dataframe['Keywords Count']:
        terms_in_doc = [term for term, c in counts.items() if c > 0]
        doc_freq.update(set(terms_in_doc))

    # Compute IDF for each term
    idf_dict = {term: math.log(N / df_t) for term, df_t in doc_freq.items()}

    return idf_dict


def calculate_tf(row):

    """Compute term frequency (TF) values for scoring keywords for a single document."""

    # Raw term counts for the document
    counts = row['Keywords Count']

    # Total word count used for normalization
    doc_len = row['Word Count']
    
    # Return TF == 0.0 if a full XML record could not be retrieved
    if doc_len == 0:
        return {k: 0.0 for k in counts}
    
    # Normalize counts into TF values
    tf_dict = {k: v / doc_len for k, v in counts.items()}

    return tf_dict