#!/usr/bin/env python3

####################################################################
# Hermes v1.1 - Open-source mining tool for open-access literature #
# 2025-08-25                                                       #
# Written by Julien Charest & Katarina Priselac                    #
####################################################################

from Bio import Entrez, Medline
from bs4 import BeautifulSoup as bs
import time
import random
from sumy.parsers.plaintext import PlaintextParser
from sumy.nlp.tokenizers import Tokenizer
from sumy.summarizers.lex_rank import LexRankSummarizer
import nltk
from nltk.tokenize import sent_tokenize
import spacy
import re
import requests

# Getting NLTK Punkt Model for Sumy
try:
    nltk.data.find("tokenizers/punkt")
except LookupError:
    nltk.download("punkt")

# Loading Spacy Model
nlp = spacy.load("en_ner_bionlp13cg_md")

# Fetch Record
def fetch_record(pmcid):
        try:
            time.sleep(random.uniform(0.3, 0.5))
            handle = Entrez.efetch(db="pmc", id=pmcid, retmode = "xml", rettype = "full")
            record = str(list(Medline.parse(handle)))
            return record
        except Exception as e:
            return None

# Summarizes Article Main Text
def summarize_text(text, sentence_count=5):
    parser = PlaintextParser.from_string(text, Tokenizer("english"))
    summarizer = LexRankSummarizer()
    summary = summarizer(parser.document, sentence_count)
    return " ".join(cleaning_summary(str(sentence)) for sentence in summary)

# Cleaning Main Text Summary
def cleaning_summary(text):
    cleaned_text = re.sub(r'\([^()]*?\)', '', text)
    cleaned_text = re.sub(r'\b(Figure|Fig\.?)\s*\d+[A-Za-z]*[^.]*', '', cleaned_text, flags=re.IGNORECASE)
    cleaned_text = re.sub(r'\s+', ' ', cleaned_text).strip()
    cleaned_text = re.sub(r'\s+([.,;:])', r'\1', cleaned_text)
    cleaned_text = re.sub(r'([.,;:])([A-Za-z])', r'\1 \2', cleaned_text)
    cleaned_text = re.sub(r';', '', cleaned_text)
    cleaned_text = re.sub(r'\.\.+', '.', cleaned_text)
    cleaned_text = re.sub(r'([a-z0-9\.,])([A-Z][a-z]+)', r'\1 \2', cleaned_text)
    cleaned_text = re.sub(r',([A-Z])', r', \1', cleaned_text)
    cleaned_text = re.sub(r'([?!])([A-Za-z])', r'\1 \2', cleaned_text)
    cleaned_text = cleaned_text[0].upper() + cleaned_text[1:]
    return cleaned_text

# Flatten NER Output List
def flatten_list(lst):
    flat = []
    for item in lst:
        if isinstance(item, list):
            flat.extend(flatten_list(item))
        else:
            flat.append(item)
    return flat

# Validate Gene with MyGene.info
def validate_gene(gene):
    q = gene.strip()

    try:
        r = requests.get("https://mygene.info/v3/query", params={"q": f"symbol:{q}", "species": "all", "size": 10, "fields": "symbol"}, timeout=8)
        r.raise_for_status()
        hits = r.json().get("hits", []) or []
    except requests.RequestException:
        return False

    q_lower = q.lower()
    for h in hits:
        sym = (h.get("symbol") or "").strip()
        if sym.lower() == q_lower:
            return True
    return False

# Validate Protein with MyGene.info
def validate_protein(protein):
    q = protein.strip()
    r = requests.get("http://mygene.info/v3/query", params={"q": f"uniprot:{q} OR symbol:{q}", "species": "all", "size": 10}, timeout=8)
    if r.status_code != 200:
        return False
    data = r.json()
    return len(data.get("hits", [])) > 0

# Validate Disease with MyDisease.info
def validate_disease(disease):
    q = disease.strip()
    r = requests.get("http://mydisease.info/v1/query", params={"q": q, "size": 10}, timeout=8)
    if r.status_code != 200:
        return False
    data = r.json()
    return len(data.get("hits", [])) > 0

# Validate Chemical with MyChemical.info
def validate_chemical(chemical):
    q = chemical.strip()
    r = requests.get("https://mychem.info/v1/query", params={"q": q, "size": 10}, timeout=8)
    if r.status_code != 200:
        return False
    data = r.json()
    return len(data.get("hits", [])) > 0

# Validate Cell with EBI OLS4
def validate_cell(cell):
    q = cell.strip()
    try:
        r = requests.get("https://www.ebi.ac.uk/ols4/api/search", params={"q": q, "ontology": "cl,clo", "rows": 25, "queryFields": "label,synonym,obo_id,short_form"}, timeout=8)
        r.raise_for_status()
        docs = r.json().get("response", {}).get("docs", []) or []
    except Exception:
        return False

    for d in docs:
        if d.get("is_obsolete") or d.get("type") != "class":
            continue
        label = (d.get("label") or "").strip().lower()
        if label == q:
            return True
        for syn in d.get("synonym") or []:
            if (syn or "").strip().lower() == q:
                return True
    return False

# Validate Organism with NCBI and EBI OLS4
def validate_organism(organism):

    _ABBR = re.compile(r"^([A-Z])\.\s*([a-z]+)$")

    q = organism.strip()

    if q.isdigit():
        r = requests.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi", params={"db":"taxonomy","term":f"{q}[TaxID]","retmode":"json","retmax":1}, timeout=8)
        return r.ok and int(r.json().get("esearchresult", {}).get("count", "0")) > 0

    qn = q.lower()
    r = requests.get("https://www.ebi.ac.uk/ols4/api/search", params={"q": q, "ontology": "ncbitaxon", "rows": 50, "queryFields": "label,synonym"}, timeout=8)
    if r.ok:
        for d in r.json().get("response", {}).get("docs", []):
            if d.get("type") != "class" or d.get("is_obsolete"):
                continue
            label = (d.get("label") or "").strip().lower()
            if label == qn:
                return True
            for syn in d.get("synonym") or []:
                if (syn or "").strip().lower() == qn:
                    return True

    m = _ABBR.match(q)
    if m:
        g_init, species = m.groups()
        r = requests.get("https://www.ebi.ac.uk/ols4/api/search", params={"q": species, "ontology": "ncbitaxon", "rows": 100, "queryFields": "label"}, timeout=8)
        if r.ok:
            for d in r.json().get("response", {}).get("docs", []):
                if d.get("type") != "class" or d.get("is_obsolete"):
                    continue
                parts = (d.get("label") or "").split()
                if len(parts) >= 2 and parts[0][:1].upper() == g_init and parts[1].lower() == species:
                    return True

    if " " not in q:
        r = requests.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi", params={"db":"taxonomy","term":f'"{q}"[Common Name]',"retmode":"json","retmax":1}, timeout=8)
        if r.ok and int(r.json().get("esearchresult", {}).get("count", "0")) > 0:
            return True
        
    return False

# Validate Tissue with EBI OLS4
def validate_tissue(tissue):
    q = tissue.strip()
    q_lower = q.lower()

    try:
        r = requests.get("https://www.ebi.ac.uk/ols4/api/search", params={"q": q, "ontology": "uberon,bto,fma", "rows": 25, "queryFields": "label,synonym,obo_id,short_form"}, timeout=8)
        r.raise_for_status()
        docs = r.json().get("response", {}).get("docs", []) or []
    except Exception:
        return False

    for d in docs:
        if d.get("is_obsolete") or d.get("type") != "class":
            continue
        label = (d.get("label") or "").strip().lower()
        if label == q_lower:
            return True
        for syn in d.get("synonym") or []:
            if (syn or "").strip().lower() == q_lower:
                return True
            
    return False

# Validate Pathway with MyPathway.info
def validate_pathway(pathway):
    q = pathway.strip()
    r = requests.get("https://mypathway.info/v1/query", params={"q": q, "size": 10}, timeout=8)
    if r.status_code != 200:
        return False
    data = r.json()
    return len(data.get("hits", [])) > 0

# Summarizing Article (Sumy) and NER (Spacy) for Article Report
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
    article = bs(str(bs_record)[introduction_pos:references_pos], "lxml")

    # Getting Figure Legends
    figs_text = []
    figs = article.find_all("fig")

    for fig in figs:
        title = fig.find("title")

        if title:
            figs_text.append(title.get_text(strip=True))
    
        else:
            ps = fig.find_all("p")

            if not ps:
                figs_text.append("NA")

            elif len(ps) == 1:
                text = ps[0].get_text(" ", strip=True)
                sentences = sent_tokenize(text)
                figs_text.append(sentences[0] if sentences else text)

            else:
                figs_text.append(ps[0].get_text(strip=True))

    # Remove Figures from Article
    for fig in article.find_all("fig"):
        fig.decompose()

    # Remove Tables from Article
    for table in article.find_all("table-wrap"):
        table.decompose()
    for table in article.find_all("table"):
        table.decompose()

    # Remove Titles from Article
    for title in article.find_all("title"):
        title.decompose()

    # Remove References from Article
    for reference in article.find_all("xref"):
        reference.decompose()
    for reference in article.find_all("ref"):
        reference.decompose()

    # Remove Supplementary Material from Article
    for suppl in article.find_all("supplementary-material"):
        suppl.decompose()

    article = article.get_text()

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

    results = [pmcid, [summary, gene_names, proteins, diseases, chemicals, cells, organisms, tissues, pathways, figs_text]]

    return results