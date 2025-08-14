#!/usr/bin/env python3
# Written by Julien Charest

from Bio import Entrez, Medline
from bs4 import BeautifulSoup as bs
import time
import random
from sumy.parsers.plaintext import PlaintextParser
from sumy.nlp.tokenizers import Tokenizer
from sumy.summarizers.lex_rank import LexRankSummarizer
import nltk
import spacy
import re
import requests
import pandas as pd

# Getting nltk model for sumy
try:
    nltk.data.find("tokenizers/punkt")
except LookupError:
    nltk.download("punkt")

# Dealing with SpaCy
nlp = spacy.load("en_ner_bionlp13cg_md")

# Fetch Record Function (already in app)
def fetch_record(pmcid):
        try:
            time.sleep(random.uniform(0.3, 0.5))
            handle = Entrez.efetch(db="pmc", id=pmcid, retmode = "xml", rettype = "full")
            record = str(list(Medline.parse(handle)))
            return record
        except Exception as e:
            return None
        
def summarize_text(text, sentence_count=5):
    parser = PlaintextParser.from_string(text, Tokenizer("english"))
    summarizer = LexRankSummarizer()
    summary = summarizer(parser.document, sentence_count)
    # Join summary sentences
    return " ".join(remove_citations_and_figures(str(sentence)) for sentence in summary)

def remove_citations_and_figures(text):
    # Remove all parentheses with citations or figure refs
    cleaned_text = re.sub(r'\([^()]*?\)', '', text)

    # Remove "Figure <number>" and related figure labels/captions appearing in-line
    cleaned_text = re.sub(r'\b(Figure|Fig\.?)\s*\d+[A-Za-z]*[^.]*', '', cleaned_text, flags=re.IGNORECASE)

    # Remove extra spaces leftover from removals
    cleaned_text = re.sub(r'\s+', ' ', cleaned_text).strip()

    # Fix spaces before punctuation (like before commas, periods)
    cleaned_text = re.sub(r'\s+([.,;:])', r'\1', cleaned_text)

    # Add a space after periods or commas if missing (glued sentences)
    cleaned_text = re.sub(r'([.,;:])([A-Za-z])', r'\1 \2', cleaned_text)

    # Remove stray semicolons leftover
    cleaned_text = re.sub(r';', '', cleaned_text)

    # Fix double periods
    cleaned_text = re.sub(r'\.\.+', '.', cleaned_text)

    # Add space between a lowercase/period/comma and an uppercase word without space (likely glued titles)
    cleaned_text = re.sub(r'([a-z0-9\.,])([A-Z][a-z]+)', r'\1 \2', cleaned_text)

    # Fix missing space after commas before uppercase letters (e.g. "NH4,Cl" -> "NH4, Cl")
    cleaned_text = re.sub(r',([A-Z])', r', \1', cleaned_text)

    # Add space after question marks or exclamation marks if glued to the next word
    cleaned_text = re.sub(r'([?!])([A-Za-z])', r'\1 \2', cleaned_text)

    # Add space after lowercase or digit followed directly by uppercase (e.g., 'identity?Below' or 'neuronThe')
    cleaned_text = re.sub(r'([a-z0-9])([A-Z])', r'\1 \2', cleaned_text)

    cleaned_text = cleaned_text[0].upper() + cleaned_text[1:]

    return cleaned_text

def flatten_list(lst):
    flat = []
    for item in lst:
        if isinstance(item, list):
            flat.extend(flatten_list(item))
        else:
            flat.append(item)
    return flat

def validate_gene(gene):
    #Validate gene name in MyGene.info
    url = "http://mygene.info/v3/query"
    params = {"q": f"symbol:{gene}", "species": "all", "size": 1}
    r = requests.get(url, params=params)
    if r.status_code != 200:
        return False
    data = r.json()
    return len(data.get("hits", [])) > 0

def validate_protein(protein):
    # Validate protein name or UniProt ID in MyGene.info
    url = "http://mygene.info/v3/query"
    params = {"q": f"uniprot:{protein} OR symbol:{protein}", "species": "all", "size": 1}
    r = requests.get(url, params=params)
    if r.status_code != 200:
        return False
    data = r.json()
    return len(data.get("hits", [])) > 0

def validate_disease(disease):
    # Validate a disease name or identifier using MyDisease.info
    url = "http://mydisease.info/v1/query"
    params = {"q": disease, "size": 1}
    r = requests.get(url, params=params)
    if r.status_code != 200:
        return False
    data = r.json()
    return len(data.get("hits", [])) > 0

def validate_chemical(chemical):
    # Validate a chemical name or identifier using MyChem.info
    url = "https://mychem.info/v1/query"
    params = {"q": chemical, "size": 1}
    r = requests.get(url, params=params)
    if r.status_code != 200:
        return False
    data = r.json()
    return len(data.get("hits", [])) > 0

def validate_cell(cell):
    # Validate a cell / cell line name or ID using EBI OLS4 (CL, CLO, EFO)
    url = "https://www.ebi.ac.uk/ols4/api/search"
    params = {"q": cell, "ontology": "cl,clo,efo", "rows": 1}
    r = requests.get(url, params=params, timeout=8)
    if r.status_code != 200:
        return False
    try:
        data = r.json()
    except ValueError:
        return False
    return (data.get("response", {}).get("numFound", 0) > 0)

def validate_organism(organism):
    # Validate an organism name or taxonomy ID using MyGene.info
    url = "https://mygene.info/v3/query"
    params = {"q": organism, "species": "all", "size": 1, "fields": "taxid,name,other_names"
    }
    r = requests.get(url, params=params)
    if r.status_code != 200:
        return False
    data = r.json()
    return len(data.get("hits", [])) > 0

def validate_tissue(tissue):
    # Validate a tissue / anatomical structure using EBI OLS4 (UBERON, BTO, FMA)
    url = "https://www.ebi.ac.uk/ols4/api/search"
    params = {"q": tissue, "ontology": "uberon,bto,fma", "rows": 1}
    r = requests.get(url, params=params, timeout=8)
    if r.status_code != 200:
        return False
    try:
        data = r.json()
    except ValueError:
        return False
    return (data.get("response", {}).get("numFound", 0) > 0)

def validate_pathway(pathway):
    # Validate a pathway name or identifier using MyPathway.info
    url = "https://mypathway.info/v1/query"
    params = {"q": pathway, "size": 1}
    r = requests.get(url, params=params)
    if r.status_code != 200:
        return False
    data = r.json()
    return len(data.get("hits", [])) > 0

def summarize_article(pmcid):
     
    record = fetch_record(pmcid)

    # Parsing Record with Beautiful Soup
    bs_record = bs(record, "lxml")
        
    # Getting Main Text
    introduction_pos = str(bs_record).find("<title>Introduction</title>")
    if introduction_pos == -1:
        introduction_pos = str(bs_record).lower().find("<title>Introduction".lower())
    references_pos = str(bs_record).find("<title>References</title>")
    if references_pos == -1:
        references_pos = str(bs_record).lower().find("<title>References".lower())
    article = bs(str(bs_record)[introduction_pos:references_pos], "lxml").get_text()

    summary = summarize_text(article, sentence_count=5)
    if len(summary) < 1:
        summary = "N/A"

    # Mentioned Biomedical Entities (SciSpacy with en_ner_bionlp13cg_md):
    doc = nlp(article)

    # - Genes ["GENE_OR_GENE_PRODUCT"]
    candidates = [ent.text for ent in doc.ents if ent.label_ in ("GENE_OR_GENE_PRODUCT")]
    flattened = set(flatten_list(candidates))
    gene_names = {gene for gene in flattened if validate_gene(gene)}
    if len(gene_names) < 1:
        gene_names = "N/A"
    else:
        gene_names = ", ".join(sorted(gene_names))

    # - Proteins ["PROTEIN"]
    candidates = [ent.text for ent in doc.ents if ent.label_ in ("PROTEIN")]
    flattened = set(flatten_list(candidates))
    proteins = {protein for protein in flattened if validate_protein(protein)}
    if len(proteins) < 1:
        proteins = "N/A"
    else:
        proteins = ", ".join(sorted(proteins))

    # - Diseases ["DISEASE"]
    candidates = [ent.text for ent in doc.ents if ent.label_ in ("DISEASE")]
    flattened = set(flatten_list(candidates))
    diseases = {disease for disease in flattened if validate_disease(disease)}
    if len(diseases) < 1:
        diseases = "N/A"
    else:
        diseases = ", ".join(sorted(diseases))

    # - Chemicals ["CHEMICAL"]
    candidates = [ent.text for ent in doc.ents if ent.label_ in ("CHEMICAL")]
    flattened = set(flatten_list(candidates))
    chemicals = {chemical for chemical in flattened if validate_chemical(chemical)}
    if len(chemicals) < 1:
        chemicals = "N/A"
    else:
        chemicals = ", ".join(sorted(chemicals))

    # - Cells ["CELL"]
    candidates = [ent.text for ent in doc.ents if ent.label_ in ("CELL")]
    flattened = set(flatten_list(candidates))
    cells = {cell for cell in flattened if validate_cell(cell)}
    if len(cells) < 1:
        cells = "N/A"
    else:
        cells = ", ".join(sorted(cells))

    # - Organisms ["ORGANISM"]
    candidates = [ent.text for ent in doc.ents if ent.label_ in ("ORGANISM")]
    flattened = set(flatten_list(candidates))
    organisms = {organism for organism in flattened if validate_organism(organism)}
    if len(organisms) < 1:
        organisms = "N/A"
    else:
        organisms = ", ".join(sorted(organisms))

    # - Tissues ["TISSUE"]
    candidates = [ent.text for ent in doc.ents if ent.label_ in ("TISSUE")]
    flattened = set(flatten_list(candidates))
    tissues = {tissue for tissue in flattened if validate_tissue(tissue)}
    if len(tissues) < 1:
        tissues = "N/A"
    else:
        tissues = ", ".join(sorted(tissues))

    # - Pathway ["PATHWAYS"]
    candidates = [ent.text for ent in doc.ents if ent.label_ in ("PATHWAYS")]
    flattened = set(flatten_list(candidates))
    pathways = {pathway for pathway in flattened if validate_pathway(pathway)}
    if len(pathways) < 1:
        pathways = "N/A"
    else:
        pathways = ", ".join(sorted(pathways))

    results = [pmcid, [summary, gene_names, proteins, diseases, chemicals, cells, organisms, tissues, pathways]]

    return results