#!/usr/bin/env python3

####################################################################
# Hermes v1.2 - Open-source mining tool for open-access literature #
# 2025-12-18                                                       #
# Written by Julien Charest & Katarina Priselac                    #
####################################################################

import os
import time
from fpdf import FPDF
import modules.adding_functions as af


def generate_pdf_report(date, terms, parsed, hits, keywords, hermes_mode, request, figures_dir, results_df, report_dir, log_file, idf_dict):

    """
    Generate and save a comprehensive HERMES PDF report summarizing the query.

    The report includes:
        - Query metadata and scoring details
        - Global statistics plots (publication years, keyword distribution, scores)
        - A summary table of all retrieved articles
        - Per-article detailed pages with metadata, abstract, summary, and entities

    Parameters
    ----------
    date : str
        Query date.
    terms : str
        Original query terms submitted to HERMES.
    parsed : int
        Number of successfully parsed articles.
    hits : int
        Total number of hits returned by the initial search.
    keywords : list[str]
        Keywords used for scoring and TF/IDF calculations.
    hermes_mode : str
        String describing HERMES scoring mode.
    request : str
        Unique identifier for the query.
    figures_dir : str
        Directory path containing pre-generated figures for the report.
    results_df : pandas DataFrame
        DataFrame with article information and scoring results.
    report_dir : str
        Directory where the final PDF report will be written.
    log_file : str
        Path to the log file where warnings/errors are appended.
    idf_dict : dict
        Mapping of keyword to IDF score used in scoring.

    """

    header_text="Hermes (v1.2) - Literature Mining Report"
    logo_path=af.ressource_path("assets/hermes_logo.png")

    # Define a custom PDF class providing a standardized layout for HERMES reports
    class PDF_Report(FPDF):
        def __init__(self):
            super().__init__()
            self.font="Arial"
            self.m=10 
            self.pw=210-2*self.m
            self.line_thickness=0
            self.height=7

        def header(self):
            """Add a standardized header to each page."""
            self.set_font(self.font, '', 10)
            self.cell(w=(self.pw/3), h=self.height, txt="", border=self.line_thickness, ln=0, align='L')
            self.cell(w=(self.pw/3), h=self.height, txt=header_text, border=self.line_thickness, ln=0, align ='C')
            self.cell(w=(self.pw/3), h=self.height, txt="{date}".format(date=date), border=self.line_thickness, ln=1, align ="R")
            
            # Attempt to draw HERMES logo
            try:
                self.image(logo_path, 10, 8, 35)
            except Exception as e:
                print("[WARNING]: Could not load logo:", logo_path, repr(e))
                with open(log_file, 'a') as f:
                    f.write("[WARNING]: Could not load logo.\n")
            self.ln(5)

        def footer(self):
            """Add a page footer with page numbering."""
            self.set_y(-15)
            self.set_font(self.font, '', 8)
            self.cell(w=(self.pw), h=self.height, txt='Page {page_number}'.format(page_number = self.page_no()), border=self.line_thickness, ln=1, align ='C')

    # -------------------------------
    # First page: title + query stats
    # -------------------------------
    pdf=PDF_Report()
    pdf.add_page()
    pdf.set_font(pdf.font, 'B', 24)
    pdf.cell(w=pdf.pw, h=25, txt="Literature Mining Report", border=pdf.line_thickness, ln=1, align='C')
    pdf.set_font(pdf.font, 'B', 12)
    
    # Query statistics block
    pdf.cell(w=pdf.pw, h=8, txt="Query Statistics", border=pdf.line_thickness, ln=1, align='L')
    pdf.set_font(pdf.font, '', 10)
    pdf.cell(w=pdf.pw, h=5, txt="Query Date: {date}".format(date=date), border=pdf.line_thickness, ln=1, align='L')
    pdf.cell(w=pdf.pw, h=5, txt='Query Terms: "{terms}"'.format(terms=terms), border=pdf.line_thickness, ln=1, align='L')
    pdf.cell(w=pdf.pw, h=5, txt='Scoring Keywords: {keywords} {mode}'.format(keywords = "; ".join(sorted(keywords)), mode=hermes_mode), border=pdf.line_thickness, ln=1, align='L')
    pdf.cell(w=pdf.pw, h=5, txt='Articles Parsed: {parsed}/{hits} ({percentage}%)'.format(parsed=int(parsed), hits=int(hits), percentage="{:.2f}".format(parsed*100/hits)), border=pdf.line_thickness, ln=1, align='L')
    pdf.cell(w=pdf.pw, h=5, txt='Keyword Rarity Score (IDF): {idf}'.format(idf = "; ".join(f"{k}: {v:.3f}" for k, v in idf_dict.items())), border=pdf.line_thickness, ln=1, align='L')
    pdf.ln(3)

    # Insert overview figures: publication years and associated keywords
    pdf.image("{figures_dir}/{request}_pubyears.png".format(figures_dir=figures_dir, request=request), w=pdf.pw)
    pdf.image("{figures_dir}/{request}_asskeywords.png".format(figures_dir=figures_dir, request=request), w=pdf.pw)
    
    # -----------------------------------------------
    # Second page: statistics summary + summary table
    # -----------------------------------------------
    pdf.add_page()
    pdf.ln(2)

    # Insert overview figures: statistics summary figure
    pdf.image("{figures_dir}/{request}_stats_summary.png".format(figures_dir=figures_dir, request=request), w=pdf.pw)

    # Summary table: Top N Results
    pdf.ln(2)
    pdf.set_font(pdf.font, 'B', 10)
    pdf.cell(w=pdf.pw, h=5, txt="Query Results", border=pdf.line_thickness, ln=1, align='C')
    pdf.ln(1)

    # Table header
    pdf.set_font(pdf.font, 'B', 9)
    pdf.table_border=1
    pdf.cell(w=(3*pdf.pw/30), h=5, txt="PMID", border=pdf.table_border, ln=0, align='C')
    pdf.cell(w=(2*pdf.pw/30), h=5, txt="Year", border=pdf.table_border, ln=0, align='C')
    pdf.cell(w=(4*pdf.pw/30), h=5, txt="Journal", border=pdf.table_border, ln=0, align='C')
    pdf.cell(w=(5*pdf.pw/30), h=5, txt="Authors", border=pdf.table_border, ln=0, align='C')
    pdf.cell(w=(11*pdf.pw/30), h=5, txt="Title", border=pdf.table_border, ln=0, align='C')
    pdf.cell(w=(3*pdf.pw/30), h=5, txt="Citations", border=pdf.table_border, ln=0, align='C')
    pdf.cell(w=(2*pdf.pw/30), h=5, txt="Score", border=pdf.table_border, ln=1, align='C')

    # Table rows: one line per article
    pdf.set_font(pdf.font, '', 8)
    for i in range(len(results_df)):
        pdf.cell(w=(3*pdf.pw/30), h=5, txt=results_df["PMID"][i], link="https://doi.org/{doi}".format(doi=results_df["DOI"][i]), border=pdf.table_border, ln=0, align='C')
        pdf.cell(w=(2*pdf.pw/30), h=5, txt=str(results_df["Year"][i]), border=pdf.table_border, ln=0, align='C')
        pdf.cell(w=(4*pdf.pw/30), h=5, txt=af.filter_string(results_df["Journal"][i][:15]), border=pdf.table_border, ln=0, align='L')
        if len(results_df["Authors"][i]) > 1:
            pdf.cell(w=(5*pdf.pw/30), h=5, txt=af.filter_string("{first_author} et al.".format(first_author=results_df["Authors"][i][0][:15])), border=pdf.table_border, ln=0, align='L')
        else:
            if len(results_df["Authors"][i]) < 1:
                pdf.cell(w=(5*pdf.pw/30), h=5, txt="N/A", border=pdf.table_border, ln=0, align='L')
            else:
                pdf.cell(w=(5*pdf.pw/30), h=5, txt=af.filter_string(results_df["Authors"][i][0][:15]), border=pdf.table_border, ln=0, align='L')
        if len(results_df["Title"][i]) > 45:
            pdf.cell(w=(11*pdf.pw/30), h=5, txt=af.filter_string("{title}...".format(title = results_df["Title"][i][:45].lower())), border=pdf.table_border, ln=0, align='L')
        else:
            pdf.cell(w=(11*pdf.pw/30), h=5, txt=af.filter_string(results_df["Title"][i]), border=pdf.table_border, ln=0, align='L')
        pdf.cell(w=(3*pdf.pw/30), h=5, txt="{citation}".format(citation=results_df["Citations"][i]), border=pdf.table_border, ln=0, align='C')
        pdf.cell(w=(2*pdf.pw/30), h=5, txt="{score}".format(score="{:.2f}".format(results_df["Score"][i])), border=pdf.table_border, ln=1, align='C')

    # ---------------------------------
    # Detailed per-article report pages
    # ---------------------------------
    counter=1
    pdf.add_page()
    for i in range(len(results_df)):
        # Start each new article on a fresh page after the first
        if counter > 1:
            pdf.add_page()
            pdf.ln(3)
            counter=1

        # Article header
        pdf.set_font(pdf.font, 'B', 10)
        pdf.cell(w=pdf.pw, h=5, txt="Query Result {i}/{len_results}".format(i=i+1, len_results=len(results_df)), border=pdf.line_thickness, align='L')
        pdf.ln(5)
        pdf.set_font(pdf.font, '', 10)
        pdf.multi_cell(w=pdf.pw, h=5, txt=af.filter_string(results_df["Title"][i]), border=pdf.line_thickness, ln=1, align='L')
        pdf.set_font(pdf.font, '', 9)
        pdf.multi_cell(w=pdf.pw, h=5, txt=af.filter_string("; ".join(results_df["Authors"][i])), border=pdf.line_thickness, ln=1, align='L')
        pdf.set_font(pdf.font, '', 8)
        pdf.cell(w=pdf.pw, h=5, txt=af.filter_string(results_df["Journal"][i]), border=pdf.line_thickness, ln=1, align='L')
        pdf.cell(w=(4*pdf.pw/30), h=5, txt="Year: {year}".format(year=results_df["Year"][i]), border=pdf.line_thickness, ln=0, align='L')
        pdf.cell(w=(4*pdf.pw/30), h=5, txt="PMID: {pmid}".format(pmid=results_df["PMID"][i]), border=pdf.line_thickness, ln=0, align='L')
        pdf.cell(w=(4*pdf.pw/30), h=5, txt="Citations: {citations}".format(citations=results_df["Citations"][i]), border=pdf.line_thickness, ln=0, align='L')
        pdf.cell(w=(4*pdf.pw/30), h=5, txt="Score: {score}".format(score="{:.2f}".format(results_df["Score"][i])), border=pdf.line_thickness, ln=0, align='L')
        pdf.cell(w=(10*pdf.pw/30), h=5, txt="DOI: {doi}".format(doi=results_df["DOI"][i]), link="https://doi.org/{doi}".format(doi=results_df["DOI"][i]), border=pdf.line_thickness, ln=1, align='L')

        # Keyword statistics (raw counts and normalized TF)
        pdf.cell(w=(4*pdf.pw/30), h=5, txt=af.filter_string("Keywords Hits: {keywords}".format(keywords=results_df["Keywords Count"][i])), border=pdf.line_thickness, ln=1, align='L')
        pdf.cell(w=(4*pdf.pw/30), h=5, txt=af.filter_string("Keywords Frequency (TF): {tf}".format(tf="; ".join(f"{k}: {v:.3f}" for k, v in results_df["TF"][i].items()))), border=pdf.line_thickness, ln=1, align='L')
        pdf.ln(2)
        
        # Abstract section
        pdf.set_font(pdf.font, 'B', 8)
        pdf.cell(w=pdf.pw, h=5, txt="Abstract:", border=pdf.line_thickness, align='L')
        pdf.ln(4)
        pdf.set_font(pdf.font, '', 8)
        pdf.multi_cell(w=pdf.pw, h=5, txt=af.filter_string(results_df["Abstract"][i]), border=pdf.line_thickness, align='L')
        pdf.ln(2)

        # Summary section
        pdf.set_font(pdf.font, 'B', 8)
        pdf.cell(w=pdf.pw, h=5, txt="Summary:", border=pdf.line_thickness, align='L')
        pdf.ln(4)
        pdf.set_font(pdf.font, '', 8)
        for j in range(len(results_df["Summary"][i][0])):
            pdf.multi_cell(w=pdf.pw, h=5, txt=af.filter_string(results_df["Summary"][i][0][j]), border=pdf.line_thickness, align='L')
            pdf.ln(1)
        pdf.ln(1)

        # Manuscript structure
        # pdf.set_font(pdf.font, 'B', 8)
        # pdf.cell(w=pdf.pw, h=5, txt="Sections:", border=pdf.line_thickness, align='L')
        # pdf.ln(4)
        # pdf.set_font(pdf.font, '', 8)
        # for j in range(len(results_df["Summary"][i][10])):
        #    pdf.multi_cell(w=pdf.pw, h=5, txt="{j}. ".format(j=j+1) + af.filter_string(results_df["Summary"][i][10][j]), border=pdf.line_thickness, align='L')
        #    pdf.ln(1)
        # pdf.ln(1)

        # Figure descriptions
        pdf.set_font(pdf.font, 'B', 8)
        pdf.cell(w=pdf.pw, h=5, txt="Figures:", border=pdf.line_thickness, align='L')
        pdf.ln(4)
        pdf.set_font(pdf.font, '', 8)
        for j in range(len(results_df["Summary"][i][9])):
            pdf.multi_cell(w=pdf.pw, h=5, txt="Figure {j}: ".format(j=j+1) + af.filter_string(results_df["Summary"][i][9][j]), border=pdf.line_thickness, align='L')
            pdf.ln(1)
        pdf.ln(1)

        # Biomedical entity lists extracted
        pdf.set_font(pdf.font, 'B', 8)
        pdf.cell(w=pdf.pw, h=5, txt="Mentioned biomedical entities:", border=pdf.line_thickness, align='L')
        pdf.ln(4)
        pdf.set_font(pdf.font, '', 8)
        pdf.multi_cell(w=pdf.pw, h=5, txt="Genes: "+af.filter_string(results_df["Summary"][i][1]), border=pdf.line_thickness, align='L')
        pdf.ln(1)
        pdf.multi_cell(w=pdf.pw, h=5, txt="Proteins: "+af.filter_string(results_df["Summary"][i][2]), border=pdf.line_thickness, align='L')
        pdf.ln(1)
        pdf.multi_cell(w=pdf.pw, h=5, txt="Diseases: "+af.filter_string(results_df["Summary"][i][3]), border=pdf.line_thickness, align='L')
        pdf.ln(1)
        pdf.multi_cell(w=pdf.pw, h=5, txt="Chemicals: "+af.filter_string(results_df["Summary"][i][4]), border=pdf.line_thickness, align='L')
        pdf.ln(1)
        pdf.multi_cell(w=pdf.pw, h=5, txt="Cells: "+af.filter_string(results_df["Summary"][i][5]), border=pdf.line_thickness, align='L')
        pdf.ln(1)
        pdf.multi_cell(w=pdf.pw, h=5, txt="Organisms: "+af.filter_string(results_df["Summary"][i][6]), border=pdf.line_thickness, align='L')
        pdf.ln(1)
        pdf.multi_cell(w=pdf.pw, h=5, txt="Tissues: "+af.filter_string(results_df["Summary"][i][7]), border=pdf.line_thickness, align='L')
        pdf.ln(1)
        pdf.multi_cell(w=pdf.pw, h=5, txt="Pathways: "+af.filter_string(results_df["Summary"][i][8]), border=pdf.line_thickness, align='L')
        counter+=1
    
    # -------------------------
    # Save PDF report to disk
    # -------------------------
    retries=2
    delay=0.5
    output_path = os.path.join(report_dir, f"{date}_{request}.pdf")
    for attempt in range(retries + 1):
        try:
            pdf.output(output_path)
        except Exception as e:
            print(f"[ERROR] PDF report generation failed on attempt {attempt+1}: {repr(e)}")
        
        time.sleep(delay)
        if os.path.exists(output_path) and os.path.getsize(output_path) > 0:
            with open(log_file, 'a') as f:
                    f.write("[STATUS] PDF report saved successfully.\n")
            return None

    print(f"[ERROR] Could not generate PDF report after {retries+1} attempts.")
    with open(log_file, 'a') as f:
                    f.write("[ERROR] Could not generate PDF report.\n")
    return None