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
- Generates reports in PDF format with summary figures
- **No installation required** – ready-to-use executables for all platforms

### Installation 

*HERMES* requires no installation. Executable for all platforms are available in [latest release](https://github.com/julien-charest/hermes/releases). Alternatively, *HERMES* can be run from source:

```r
# Clone the repository:
git clone https://github.com/julien-charest/hermes.git
cd hermes

# Install dependencies manually:
pip install biopython pandas beautifulsoup4 fpdf lxml

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