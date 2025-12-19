#!/usr/bin/env python3

####################################################################
# Hermes v1.1 - Open-source mining tool for open-access literature #
# 2025-12-18                                                       #
# Written by Julien Charest & Katarina Priselac                    #
####################################################################

from bs4 import BeautifulSoup as bs
import re
import time
import random
import requests
from Bio import Entrez
from transformers import pipeline, T5Tokenizer, T5ForConditionalGeneration
import spacy
from transformers.utils import logging as hf_logging
import warnings


##################################
# Configure Performance Settings #
##################################

# LLM model configuration
#   Select the local or remote Hugging Face-supported model used for summarization.
LLM_MODEL_NAME = "google/flan-t5-large"

# NER model configuration
#   Select the local or remote Spacy-supported model for named entity recognition.
NER_MODEL_NAME = "en_ner_bionlp13cg_md"

# Compute device configuration
## device = -1               # CPU
## device = 0                # First CUDA GPU
device = 0


################################
# Loading Summarization Models #
################################

# Initialize the Hugging Face T5-based LLM summarization pipeline and suppress non-critical warnings
print("[STATUS] Loading {MODEL_NAME} model for LLM summarization...".format(MODEL_NAME = LLM_MODEL_NAME))
hf_logging.set_verbosity_error()
_tokenizer = T5Tokenizer.from_pretrained(LLM_MODEL_NAME, legacy = False)
_model = T5ForConditionalGeneration.from_pretrained(LLM_MODEL_NAME)
_SUMMARIZER = pipeline("summarization", model = _model, tokenizer = _tokenizer, device = device)

# Initialize the spaCy-based NER pipeline and suppress known non-critical warnings
print("[STATUS] Loading {MODEL_NAME} model for NER...".format(MODEL_NAME = NER_MODEL_NAME))
warnings.filterwarnings("ignore", message=r"Possible set union at position .*", category=FutureWarning)
nlp = spacy.load("en_ner_bionlp13cg_md")


###############################
# HERMES Summarization Module #
###############################

def summarize_sections(sections, min_new_tokens = 100, max_new_tokens = 2500, batch_size = 5):

    """Summarize each section using a prompt-aware batch call."""

    instruction = ("Summarize the following scientific text in 3–6 sentences. Only include information explicitly mentioned. Maintain scientific accuracy and correct formatting.\nTEXT:\n")
    titles = []
    inputs = []

    for title, text in sections.items():
        if text and text.strip():
            titles.append(title)

            # Build prompt for each section
            cleaned = clean_text(text.strip())
            formatted_input = f"{instruction}{cleaned}"
            inputs.append(formatted_input)

    # Nothing to summarize
    if not inputs:
        return {title: "" for title in sections}

    # Run model ONCE for the batch
    results = _SUMMARIZER(
        inputs,
        truncation=True,
        min_new_tokens=min_new_tokens,
        max_new_tokens=max_new_tokens,
        num_beams=1,
        no_repeat_ngram_size=3,
        do_sample=False,
        batch_size=batch_size,
    )

    summaries = {
        title: result.get("summary_text", "").strip()
        for title, result in zip(titles, results)
    }

    # Guarantee all original keys exist
    for title in sections.keys():
        summaries.setdefault(title, "")

    return summaries


def extract_sections(bs_record):

    """
    Split an XML scientific article into structured sections based on headings.
    Returns a dict like:
      {"Introduction": "text...", "Methods": "text...", ...}
    """

    sections = {}

    # Find all <title> tags
    title_tags = bs_record.find_all("title")

    # If no <title> exists, fallback to full text
    if not title_tags:
        return {"Untitled": bs_record.get_text(" ", strip=True)}

    # Loop through each <title>
    for i, title_tag in enumerate(title_tags):
        section_name = title_tag.get_text(strip=True)

        # Determine next title to know when section ends
        next_title = title_tags[i+1] if i+1 < len(title_tags) else None

        content = []
        skip_first = True  # ensures we skip the title text itself

        # Iterate through elements after this title until the next one
        for elem in title_tag.next_elements:

            # Stop if we reached the next section title
            if elem == next_title:
                break

            # Skip the first text node (belongs to <title>)
            if skip_first:
                skip_first = False
                continue

            # Extract text content
            if hasattr(elem, "get_text"):
                text = elem.get_text(strip=True)
                if text and text != section_name:   # safety check
                    content.append(text)

        sections[section_name] = " ".join(content).strip()

    return sections


# Summarizing Article (LLM) and NER (Spacy) for Article Report
def summarize_article(pmcid):

    """
    Generate an article-level report by combining LLM-based section summarization
    with biomedical named entity recognition (NER).

    This function retrieves a full-text PMC article, preprocesses and filters its
    content, produces section-wise summaries using a Hugging Face T5 model, and
    extracts key biomedical entities (genes, proteins, diseases, chemicals, etc.)
    using a spaCy/SciSpacy NER pipeline. The resulting summary and extracted entities
    are returned in a structured format for inclusion in the final report.
    """

    # Retrieve the PMC article record in XML format
    record = fetch_record_xml(pmcid)
    if record is None:
         raise RuntimeError(f"Entrez fetch failed for {pmcid}")
    
    # Parse the XML record into a BeautifulSoup object for structured content access
    try:
        article = bs(record, "lxml-xml")
    except Exception as e:
        raise RuntimeError(f"Could not parse XML for {pmcid}: {e}")

    # Extract figure legends separately for inclusion in the final report
    figs_text = get_figures(article)

    # Remove table content to avoid summarizing dense tabular data
    for table in article.find_all("table-wrap"):
        table.decompose()
    for table in article.find_all("table"):
        table.decompose()

    # Remove reference markers and reference sections to reduce citation noise
    for reference in article.find_all("xref"):
        reference.decompose()
    for reference in article.find_all("ref"):
        reference.decompose()

    # Remove supplementary material to focus on the main manuscript content
    for suppl in article.find_all("supplementary-material"):
        suppl.decompose()

    # Normalize inline formatting tags to plain text for cleaner downstream processing
    INLINE_TAGS = ["italic", "bold", "sup", "sub", "xref", "ext-link"]
    for tag in article.find_all(INLINE_TAGS):
        tag.unwrap()

    # Extract article sections and retain section titles for reporting
    sections = extract_sections(article)
    sections_title = list(sections.keys())
    sections_title = [s for s in sections_title if s not in figs_text]

    # Filter to core narrative sections typically used for scientific summaries
    allowed_keys = {"introduction", "background", "results", "discussion", "results and discussion", "conclusion", "conclusions",
                    "concluding remarks", "perspectives", "conclusion and perspectives", "reporting summary"}
    filtered_sections = {k: v for k, v in sections.items() if k.lower() in allowed_keys}
    
    # Generate section-level summaries with the LLM and combine into a single report block
    summaries = summarize_sections(filtered_sections)
    summary = dict_to_block_text(summaries)
    
    # Run biomedical NER over the cleaned full text to extract mentioned entities
    doc = nlp(article.get_text())

    # Extract and validate gene mentions from NER output
    candidates = [ent.text for ent in doc.ents if ent.label_ in ("GENE_OR_GENE_PRODUCT")]
    flattened = set(flatten_list(candidates))
    gene_names = {gene for gene in flattened if validate_gene(gene)}
    if len(gene_names) < 1:
        gene_names = "N/A"
    else:
        gene_names = ", ".join(sorted(gene_names))

    # Extract and validate protein mentions from NER output
    candidates = [ent.text for ent in doc.ents if ent.label_ in ("PROTEIN")]
    flattened = set(flatten_list(candidates))
    proteins = {protein for protein in flattened if validate_protein(protein)}
    if len(proteins) < 1:
        proteins = "N/A"
    else:
        proteins = ", ".join(sorted(proteins))

    # Extract and validate disease mentions from NER output
    candidates = [ent.text for ent in doc.ents if ent.label_ in ("DISEASE")]
    flattened = set(flatten_list(candidates))
    diseases = {disease for disease in flattened if validate_disease(disease)}
    if len(diseases) < 1:
        diseases = "N/A"
    else:
        diseases = ", ".join(sorted(diseases))

    # Extract and validate chemical mentions from NER output
    candidates = [ent.text for ent in doc.ents if ent.label_ in ("CHEMICAL")]
    flattened = set(flatten_list(candidates))
    chemicals = {chemical for chemical in flattened if validate_chemical(chemical)}
    if len(chemicals) < 1:
        chemicals = "N/A"
    else:
        chemicals = ", ".join(sorted(chemicals))

    # Extract and validate cell mentions from NER output
    candidates = [ent.text for ent in doc.ents if ent.label_ in ("CELL")]
    flattened = set(flatten_list(candidates))
    cells = {cell for cell in flattened if validate_cell(cell)}
    if len(cells) < 1:
        cells = "N/A"
    else:
        cells = ", ".join(sorted(cells))

    # Extract and validate organism mentions from NER output
    candidates = [ent.text for ent in doc.ents if ent.label_ in ("ORGANISM")]
    flattened = set(flatten_list(candidates))
    organisms = {organism for organism in flattened if validate_organism(organism)}
    if len(organisms) < 1:
        organisms = "N/A"
    else:
        organisms = ", ".join(sorted(organisms))

    # Extract and validate tissue mentions from NER output
    candidates = [ent.text for ent in doc.ents if ent.label_ in ("TISSUE")]
    flattened = set(flatten_list(candidates))
    tissues = {tissue for tissue in flattened if validate_tissue(tissue)}
    if len(tissues) < 1:
        tissues = "N/A"
    else:
        tissues = ", ".join(sorted(tissues))

    # Extract and validate pathway mentions from NER output
    candidates = [ent.text for ent in doc.ents if ent.label_ in ("PATHWAYS")]
    flattened = set(flatten_list(candidates))
    pathways = {pathway for pathway in flattened if validate_pathway(pathway)}
    if len(pathways) < 1:
        pathways = "N/A"
    else:
        pathways = ", ".join(sorted(pathways))

    # Assemble the summary, extracted entities, and auxiliary metadata for downstream reporting
    results = [summary, gene_names, proteins, diseases, chemicals, cells, organisms, tissues, pathways, figs_text, sections_title]
    return results


def dict_to_block_text(dict):

    """Convert a dictionary of section summaries into labeled text blocks for report assembly."""

    return [f"{key}: {value}" for key, value in dict.items()]


def clean_text(text):

    """
    Normalize and sanitize free-text content by removing formatting artifacts,
    noisy punctuation, and unsupported characters while preserving scientific notation.
    """

    # Return empty string if input text is missing
    if not text:
        return ""

    # Normalize whitespace by replacing newlines and tabs with spaces
    text = text.replace("\n", " ").replace("\t", " ")

    # Remove duplicated punctuation sequences
    text = re.sub(r"[;:/|><]{2,}", " ", text)

    # Remove repeated punctuation blocks  
    text = re.sub(r"[\.\,\;\:\-\(\)\'\"]{3,}", " ", text)

    # Remove isolated punctuation characters surrounded by whitespace
    text = re.sub(r"\s[^\w\s]\s", " ", text)

    # Collapse multiple consecutive spaces into a single space
    text = re.sub(r"\s{2,}", " ", text)

    # Remove unsupported characters while preserving common scientific symbols
    text = re.sub(r"[^a-zA-Z0-9µα-ωΑ-Ω\-\.,;:\(\)\/ ]", "", text)

    # Trim leading and trailing whitespace
    return text.strip()


def flatten_list(lst):

    """Recursively flatten a nested list structure into a single-level list."""

    # Initialize container for flattened elements
    flat = []

    # Recursively expand nested lists
    for item in lst:
        if isinstance(item, list):
            flat.extend(flatten_list(item))
        else:
            flat.append(item)
    return flat


def validate_gene(gene):

    """Validate a gene symbol by querying the MyGene.info API and checking for an exact symbol match."""

    # Normalize and sanitize the candidate gene symbol
    q = gene.strip()
    q_lower = q.lower()

    # Query MyGene.info for matching gene symbols across species
    try:
        r = requests.get("https://mygene.info/v3/query", params={"q": f"symbol:{q}", "species": "all", "size": 10, "fields": "symbol"}, timeout=8)
        r.raise_for_status()
        hits = r.json().get("hits", []) or []
    except requests.RequestException:
        return False

    # Perform case-insensitive exact matching against returned gene symbols
    for h in hits:
        sym = (h.get("symbol") or "").strip()
        if sym.lower() == q_lower:
            return True
    return False


def validate_protein(protein):

    """
    Validate a protein identifier by querying the MyGene.info API and checking
    for an exact case-insensitive match to UniProt accessions or symbols.
    """

    # Normalize and sanitize the candidate protein identifier
    q = protein.strip()
    q_lower = q.lower()

    # Query MyGene.info for matching UniProt accessions or protein symbols
    try:
        r = requests.get("https://mygene.info/v3/query", params={"q": f"uniprot:{q} OR symbol:{q}", "species": "all", "size": 10}, timeout=8)
        r.raise_for_status()
        hits = r.json().get("hits", []) or []
    except requests.RequestException:
        return False

    # Perform case-insensitive exact matching against returned identifiers
    for h in hits:
        symbol = (h.get("symbol") or "").strip().lower()
        uniprot = (h.get("uniprot", {}).get("Swiss-Prot") or "")
        if isinstance(uniprot, str):
            uniprot = [uniprot]
        uniprot = [u.lower() for u in uniprot]

        if q_lower == symbol or q_lower in uniprot:
            return True

    return False


def validate_disease(disease):

    """
    Validate a disease name by querying the MyDisease.info API and checking
    for a case-insensitive exact match against returned disease labels.
    """

    # Normalize and sanitize the candidate disease name
    q = disease.strip()
    q_lower = q.lower()

    # Query MyDisease.info for matching disease entries
    try:
        r = requests.get("https://mydisease.info/v1/query", params={"q": q, "size": 10}, timeout=8)
        r.raise_for_status()
        hits = r.json().get("hits", []) or []
    except requests.RequestException:
        return False

    # Perform case-insensitive exact matching against disease names and synonyms
    for h in hits:
        name = (h.get("name") or "").strip().lower()
        if name == q_lower:
            return True

        # Also check known synonyms, if present
        synonyms = h.get("synonyms") or []
        if any(q_lower == s.strip().lower() for s in synonyms):
            return True
        
    return False


def validate_chemical(chemical):

    """
    Validate a chemical name by querying the MyChem.info API and checking
    for a case-insensitive exact match against returned chemical names or synonyms.
    """

    # Normalize and sanitize the candidate chemical name
    q = chemical.strip()
    q_lower = q.lower()

    try:
        # Query MyChem.info for matching chemical entries
        r = requests.get("https://mychem.info/v1/query", params={"q": q, "size": 10}, timeout=8)
        r.raise_for_status()
        hits = r.json().get("hits", []) or []
    except requests.RequestException:
        return False

    # Perform case-insensitive exact matching against chemical names and synonyms when available
    for h in hits:
        name = (h.get("name") or "").strip().lower()
        if name == q_lower:
            return True

        synonyms = h.get("synonyms") or []
        if any(q_lower == s.strip().lower() for s in synonyms):
            return True

    return False



def validate_cell(cell):

    """
    Validate a cell type name by querying the EBI OLS4 API and checking for a
    case-insensitive exact match against ontology labels or synonyms (CL/CLO).
    """

    # Normalize and sanitize the candidate cell type label
    q = cell.strip()
    q_lower = q.lower()

    # Query OLS4 for matching Cell Ontology (CL) and Cell Line Ontology (CLO) terms
    try:
        r = requests.get("https://www.ebi.ac.uk/ols4/api/search", params={"q": q, "ontology": "cl,clo", "rows": 25, "queryFields": "label,synonym,obo_id,short_form"}, timeout=8)
        r.raise_for_status()
        docs = r.json().get("response", {}).get("docs", []) or []
    except Exception:
        return False

    # Scan candidate ontology terms and match against preferred labels and synonyms
    for d in docs:
        # Skip obsolete entries and non-class results
        if d.get("is_obsolete") or d.get("type") != "class":
            continue

        # Validate against the preferred term label
        label = (d.get("label") or "").strip().lower()
        if label == q_lower:
            return True

        # Validate against any listed synonyms
        for syn in d.get("synonym") or []:
            if (syn or "").strip().lower() == q_lower:
                return True

    return False



def validate_organism(organism):

    """
    Validate an organism name using NCBI Taxonomy and EBI OLS4 (NCBITaxon),
    supporting TaxID input, exact label/synonym matching, and abbreviated binomials.
    """
    
    # Pattern to detect abbreviated binomials (e.g., "T. reesei")
    _ABBR = re.compile(r"^([A-Z])\.\s*([a-z]+)$")

    # Normalize and sanitize the candidate organism string
    q = organism.strip()

    # Validate numeric inputs as NCBI Taxonomy IDs (TaxID)
    if q.isdigit():
        r = requests.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi", params={"db": "taxonomy", "term": f"{q}[TaxID]", "retmode": "json", "retmax": 1}, timeout=8)
        return r.ok and int(r.json().get("esearchresult", {}).get("count", "0")) > 0

    # Attempt exact case-insensitive matching against NCBITaxon labels and synonyms via OLS4
    qn = q.lower()
    r = requests.get("https://www.ebi.ac.uk/ols4/api/search", params={"q": q, "ontology": "ncbitaxon", "rows": 50, "queryFields": "label,synonym"}, timeout=8)
    if r.ok:
        for d in r.json().get("response", {}).get("docs", []):
            # Skip obsolete entries and non-class results
            if d.get("type") != "class" or d.get("is_obsolete"):
                continue

            # Validate against the preferred label
            label = (d.get("label") or "").strip().lower()
            if label == qn:
                return True

            # Validate against any listed synonyms
            for syn in d.get("synonym") or []:
                if (syn or "").strip().lower() == qn:
                    return True

    # Handle abbreviated genus format (e.g., "T. reesei") by matching genus initial + species epithet
    m = _ABBR.match(q)
    if m:
        g_init, species = m.groups()
        r = requests.get("https://www.ebi.ac.uk/ols4/api/search", params={"q": species, "ontology": "ncbitaxon", "rows": 100, "queryFields": "label"}, timeout=8)
        if r.ok:
            for d in r.json().get("response", {}).get("docs", []):
                # Skip obsolete entries and non-class results
                if d.get("type") != "class" or d.get("is_obsolete"):
                    continue

                # Validate abbreviated binomial by checking genus initial and species name
                parts = (d.get("label") or "").split()
                if len(parts) >= 2 and parts[0][:1].upper() == g_init and parts[1].lower() == species:
                    return True

    # If no binomial is provided, try matching the term as a common name in NCBI Taxonomy
    if " " not in q:
        r = requests.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi", params={"db": "taxonomy", "term": f"\"{q}\"[Common Name]", "retmode": "json", "retmax": 1}, timeout=8)
        if r.ok and int(r.json().get("esearchresult", {}).get("count", "0")) > 0:
            return True

    return False


def validate_tissue(tissue):

    """
    Validate a tissue name by querying the EBI OLS4 API and checking for a
    case-insensitive exact match against ontology labels or synonyms
    (UBERON, BTO, FMA).
    """

    # Normalize and sanitize the candidate tissue name
    q = tissue.strip()
    q_lower = q.lower()

    # Query OLS4 for matching tissue-related ontology terms
    try:
        r = requests.get("https://www.ebi.ac.uk/ols4/api/search", params={"q": q, "ontology": "uberon,bto,fma", "rows": 25, "queryFields": "label,synonym,obo_id,short_form"}, timeout=8)
        r.raise_for_status()
        docs = r.json().get("response", {}).get("docs", []) or []
    except Exception:
        return False

    # Scan candidate ontology terms and validate against labels and synonyms
    for d in docs:
        # Skip obsolete entries and non-class results
        if d.get("is_obsolete") or d.get("type") != "class":
            continue

        # Validate against the preferred term label
        label = (d.get("label") or "").strip().lower()
        if label == q_lower:
            return True

        # Validate against any listed synonyms
        for syn in d.get("synonym") or []:
            if (syn or "").strip().lower() == q_lower:
                return True

    return False


def validate_pathway(pathway):

    """
    Validate a biological pathway name by querying the MyPathway.info API and
    checking for a case-insensitive exact match against returned pathway names
    or synonyms.
    """

    # Normalize and sanitize the candidate pathway name
    q = pathway.strip()
    q_lower = q.lower()

    # Query MyPathway.info for matching pathway entries
    try:
        r = requests.get("https://mypathway.info/v1/query", params={"q": q, "size": 10}, timeout=8)
        r.raise_for_status()
        hits = r.json().get("hits", []) or []
    except requests.RequestException:
        return False

    # Perform case-insensitive exact matching against pathway names and synonyms
    for h in hits:
        name = (h.get("name") or "").strip().lower()
        if name == q_lower:
            return True

        # Check pathway synonyms when available
        synonyms = h.get("synonyms") or []
        if any(q_lower == (s or "").strip().lower() for s in synonyms):
            return True

    return False


def get_figures(article):

    """
    Extract figure titles or captions from a PMC article for inclusion in the report.
    """

    # Initialize container for extracted figure text
    figs_text = []

    # Locate all figure elements in the article
    figs = article.find_all("fig")

    # Iterate over each figure to extract caption information
    for fig in figs:
        # Retrieve the caption block associated with the figure
        caption = fig.find("caption")

        # Handle figures without captions
        if not caption:
            figs_text.append("N/A")
            continue

        # Prefer explicit figure titles when available
        title = caption.find("title")
        if title:
            figs_text.append(title.get_text(" ", strip=True))
            continue

        # Fallback to paragraph-based captions
        ps = caption.find_all("p")

        # Handle empty caption paragraphs
        if not ps:
            figs_text.append("N/A")

        # Use the single paragraph if only one is present
        elif len(ps) == 1:
            text = ps[0].get_text(" ", strip=True)
            figs_text.append(text)

        # Otherwise, use the first paragraph as a concise caption
        else:
            figs_text.append(ps[0].get_text(" ", strip=True))

    # Remove placeholder entries when at least one valid caption was extracted
    if len(figs_text) > 1:
        figs_text = [t for t in figs_text if t != "N/A"]

    return figs_text


def fetch_record_xml(pmcid):
    
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
        record = handle.read()
        handle.close()

        return record

    except Exception:
        return None