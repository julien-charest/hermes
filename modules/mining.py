#!/usr/bin/env python3

####################################################################
# Hermes v1.2 - Open-source mining tool for open-access literature #
# 2025-12-18                                                       #
# Written by Julien Charest & Katarina Priselac                    #
####################################################################

import time
import random
from Bio import Entrez, Medline
from bs4 import BeautifulSoup as bs
from urllib.request import urlopen

def fetch_record(pmcid):
    
    """Retrieve a full XML record from PubMed Central for a given PMCID."""

    # Determine request pacing depending on whether an API key is present
    try:
        if Entrez.api_key:
            # API key available: higher rate allowed (~10 req/sec)
            delay = random.uniform(0.08, 0.15)
        else:
            # No API key: default NCBI rate limit (~3 req/sec)
            delay = random.uniform(0.30, 0.50)

        time.sleep(delay)

        # Retrieve the full PMC XML metadata record
        handle = Entrez.efetch(db="pmc", id=pmcid, retmode="xml", rettype="full")
        record = str(list(Medline.parse(handle)))

        return record

    except Exception:
        return None


def mine_pmcid(pmcid, keywords):

    """
    Retrieve and parse a PMC article, extracting metadata, text, and keyword stats.

    This function:
        - Fetches the full XML record from PMC for a given PMCID
        - Parses core metadata (title, journal, year, authors, DOI, keywords)
        - Extracts abstract and main text
        - Computes word counts and keyword occurrence statistics
        - Retrieves citation counts via PubMed

    Parameters
    ----------
    pmcid : str
        PubMed Central ID.
    keywords : list[str]
        List of scoring keywords used for counting occurrences in the article text.

    """

    # Fetch full PMC XML record
    record = fetch_record(pmcid)
    if record is None:
         raise RuntimeError(f"Entrez fetch failed for {pmcid}")
        
    # Parse core metadata from the XML using BeautifulSoup
    try:    
        bs_record = bs(record, "lxml")
        article_title = bs_record.find("article-title").get_text()
        journal_title = bs_record.find("journal-title").get_text()
        associated_keywords = [i.get_text().lower() for i in bs_record.find_all("kwd")]
        pub_year = bs_record.find("pub-date").find("year").get_text()
    
    except Exception as e:
        raise RuntimeError(f"Could not parse XML for {pmcid}: {e}")

    # Extract abstract text (fallback to "NA" if not present or malformed)
    try:
        abstract = bs_record.find("abstract").find("p").get_text()

    except AttributeError:
        abstract = "NA"

    except NameError:
        abstract = "NA"

    except:
        if len(bs_record.find_all("abstract")) > 0:
            abstract = "NA"
            for i in bs_record.find_all("abstract"):
                if len(i.get_text()) > len(abstract):
                    abstract = i
                    abstract = abstract.find("p").get_text()
        else:
            abstract = "NA"

    # Extract PMID and DOI, if available
    article_ids = bs_record.find_all('article-id')
    pmid = None
    doi = None
    for id in article_ids:
        if "pmid" in str(id):
            pmid = id.get_text()
        if "doi" in str(id):
            doi = id.get_text()

    # Extract author list
    contribs = bs_record.find_all("contrib")
    authors = []
    try:
        for author in contribs:
            authors.append("{surname}, {name}".format(surname = author.find("name").find("surname").get_text(), name = author.find("name").find("given-names").get_text()))
    except:
        authors = ["NA"]

    # Retrieve number of PubMed citations
    n_citations = get_citations(pmid)

    # Identify main text region
    introduction_pos = str(bs_record).find("<title>Introduction</title>")
    if introduction_pos == -1:
        introduction_pos = str(bs_record).lower().find("<title>Introduction".lower())
    references_pos = str(bs_record).find("<title>References</title>")
    if references_pos == -1:
        references_pos = str(bs_record).lower().find("<title>References".lower())
    article = bs(str(bs_record)[introduction_pos:references_pos], "lxml").get_text().lower()

    # Count total number of words in the main text
    word_count = len(article.strip().split())
    
    # Count keyword occurrences in the main text (preferred), fallback to abstract when needed
    if word_count > 0:
        keywords_count = keyword_counter(keywords, article)
    else:
        keywords_count = keyword_counter(keywords, abstract)

    # Return article mining results
    return {"PMCID": pmcid,
            "PMID": pmid,
            "Title": article_title,
            "Year": int(pub_year),
            "Journal": journal_title,
            "DOI": doi,
            "Authors": authors,
            "Associated Keywords": associated_keywords,
            "Citations": n_citations,
            "Keywords Count": keywords_count,
            "Score": 0,
            "TF": None,
            "Word Count": word_count,
            "Abstract": abstract}


def get_citations(pmid):
    
    """
    Retrieve the number of PubMed citations ("Cited by") for a given PubMed ID using the Entrez E-utilities elink endpoint.

    Parameters
    ----------
    pmid : str or int
        PubMed ID of the article to query.

   """

    # Normalize ID as a clean string
    pmid = str(pmid).strip()
    n_citations = 0

    # Retry loop to handle transient network or API errors
    for attempt in range(1, 6):
        try:
            # Determine request pacing depending on whether an API key is present
            if Entrez.api_key:
                # API key available: higher rate allowed (~10 req/sec)
                delay = random.uniform(0.08, 0.15)
            else:
                # No API key: default NCBI rate limit (~3 req/sec)
                delay = random.uniform(0.30, 0.50)

            time.sleep(delay)

            # Query PubMed citation links (pubmed_pubmed_citedin)
            handle = Entrez.elink(dbfrom="pubmed", db="pubmed", linkname="pubmed_pubmed_citedin", id=pmid)
            record = Entrez.read(handle)
            handle.close()

            # Extract citing PMIDs if available
            linksetdb = record[0].get("LinkSetDb", [])
            if linksetdb:
                links = linksetdb[0].get("Link", [])
                n_citations = len(links)
            else:
                n_citations = 0

            break

        except Exception as e:
            print(f"[DEBUG] Unexpected error for PMID {pmid}: {e}")
            break
    
    return n_citations


def keyword_counter(keywords, article):

    """
    Count occurrences of a list of keywords within a text article.

    Performs a case-insensitive search and returns a dictionary mapping each keyword to the number of times it appears in the article.

    Parameters
    ----------
    keywords : list[str]
        List of keywords to search for.
    article : str
        Full article text in which keyword occurrences are counted.

  """

    # Initialize dictionary
    keyword_count = {}

    # Normalize text once for case-insensitive matching
    text = article.lower()

    # Count occurrences for each keyword
    for keyword in keywords:
        keyword_count[keyword] = text.count(keyword.lower())

    return keyword_count