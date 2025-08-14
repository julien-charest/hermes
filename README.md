![Hermes Logo](assets/hermes_logo.png "Hermes Logo")

*HERMES* is a python-based, open-source literature mining tool designed to parse full-text, open-access articles from PubMed Central (PMC) and rank them based on:

- Occurrence of user-defined keywords  
- Recency of publication  
- Number of citations  

Its goal is to help researchers efficiently identify relevant scientific literature across vast PMC archives.

Authors: Julien Charest & Katarina Priselac

Key features:

- Parses complete articles from PubMed Central (Open Access subset)
- Ranks results using a customizable scoring algorithm
- Clean and intuitive GUI (built with Tkinter)
- Summarization of articles using Sumy with LexRank algorithm
- Biomedical entity extraction (genes, proteins, diseases, chemicals, cells, organisms, tissues, pathways) using SciSpaCy and *en_ner_bionlp13cg_md* model
- Generates reports in PDF format with summary figures

### Installation 

To install *HERMES*:

```r
# Clone the repository:
git clone https://github.com/julien-charest/hermes.git
cd hermes

# Install dependencies manually:
pip install biopython pandas beautifulsoup4 fpdf lxml matplotlib sumy nltk spacy scispacy requests

# Download scispacy model:
pip install https://s3.amazonaws.com/allenai-scispacy/models/en_ner_bionlp13cg_md-0.5.0.tar.gz

# Run the application:
python ./hermes.py
```

### Using *HERMES*

To use *HERMES* and generate a PDF report with ranked articles for your query:

- Specify a request title
- Enter your email address (for access to Entrez server)
- Define your Pubmed query (Pubmed filters are supported) (e.g. "ASE asymmetry in C. elegans")
- Enter scoring keywords (e.g. "ASE, lsy-6, che-1")
- Specify if strict mode is desired (all keywords must be included in the main text; powerful to discard article in which keywords are not mentioned together)
- Enter a number of results to be included in final report
- Submit query

![Hermes App](assets/app_screenshot.png "Hermes App")