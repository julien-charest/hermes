#!/usr/bin/env python3

####################################################################
# Hermes v1.1 - Open-source mining tool for open-access literature #
# 2025-08-25                                                       #
# Written by Julien Charest & Katarina Priselac                    #
####################################################################

from os import path, getcwd
import sys
import ssl
ssl._create_default_https_context = ssl._create_unverified_context
import time
import tkinter as tk
from tkinter import ttk, filedialog
import tkinter.font as tkFont
from datetime import datetime
from Bio import Entrez, Medline
import pandas as pd
import math
from bs4 import BeautifulSoup as bs
import add_functions as af
import summarize as summarize
from fpdf import FPDF
import numpy as np
from concurrent.futures import ThreadPoolExecutor, as_completed
import random

#####################
# Taking User Input #
#####################

# Request Title
def title_request_button():
    request_title.set(str(title_entry.get()))
    if len(request_title.get()) > 0:
        request_title_ok.set(True)
        title_entry_field["text"] = request_title.get()
    else:
        request_title_ok.set(False)
        title_entry_field["text"] = "Invalid input. Please define a request title."

# Email (mandatory for PMC requests via Entrez)
def email_button():
    email.set(str(email_entry.get()))
    if all(x in email.get() for x in ["@", "."]):
        email_ok.set(True)
        email_entry_field["text"] = email.get()
    else:
        email_ok.set(False)
        email_entry_field["text"] = "Invalid input. Please enter a valid email address: "

# Terms for PMC Query (eg. "cellulase production in Trichoderma reesei")
def terms_button():
    terms.set(str(terms_entry.get()))
    if type(terms.get()) == str and len(terms.get()) > 0:
        terms_ok.set(True)
        terms_entry_field["text"] = terms.get()
    else:
        terms_ok.set(False)
        terms_entry_field["text"] = "Invalid input. Please define PubMed PMC query terms."

# Parsing Keywords (eg. ["ace3", "xyr1"])
def keywords_button():
    keywords = []
    keywords_temp2 = keywords_entry.get().strip("'")
    if len(keywords_temp2) > 0:
        keywords_temp3 = keywords_temp2.split(",")
        for i in keywords_temp3:
            keywords.append(i.strip(" "))
        keywords_ok.set(True)
        if len(keywords) == 1:
            keywords_temp.set(keywords[0])
            keywords_entry_field["text"] = keywords_temp.get()
        else:
            keywords_list_text = ""
            for j in range(len(keywords) - 1):
                keywords_list_text += "{keyword}, ".format(keyword = keywords[j])
            keywords_list_text += keywords[-1]
            keywords_temp.set(keywords_list_text)
            keywords_entry_field["text"] = keywords_temp.get()
    else:
        keywords_temp.set("")
        keywords_ok.set(False)
        keywords_entry_field["text"] = "Invalid input. Please define a valid keywords list."

# Strict Mode (All Keywords Required)
def mandatory_keywords_check():
    if keywords_check.get() == 1:
        keywords_mandatory.set(True)
    else:
        keywords_mandatory.set(False)

# Number of Results to Include in Report
def nresults_button():
    try:
        if int(nresults_entry.get()) > 0:
            nresults.set(int(nresults_entry.get()))
            nresults_entry_field["text"] = str(nresults.get())
            nresults_ok.set(True)
        else:
            nresults_ok.set(False)
            nresults_entry_field["text"] = "Invalid input. Please define a valid number of articles."
    except ValueError:
        nresults_ok.set(False)
        nresults_entry_field["text"] = "Invalid input. Please define a valid number of articles."

# Validating User Input
def submit_query():
    if request_title_ok.get() and email_ok.get() and terms_ok.get() and keywords_ok.get() and nresults_ok.get():
        submit_button.config(state="disabled")
        submit_label["text"] = "Please wait..."
        progress_label0["text"] = ""
        progress_label["text"] = ""
        progress_label2["text"] = ""
        my_dir.set(filedialog.askdirectory())
        try:
            launch_script()
        finally:
            submit_button.config(state="normal")
            window.destroy()
            sys.exit(0)
    else:
        submit_label["text"] = "Invalid query. Please revise query arguments."

# Launch Literature Miner

def launch_script():
    request_var = request_title.get()
    email_var = email.get()
    terms_var = terms.get()
    keywords_temp_var = keywords_temp.get()
    mandatory_var = keywords_mandatory.get()
    nresults_var = nresults.get()
    literature_miner(request_var, email_var, terms_var, keywords_temp_var, mandatory_var, nresults_var)

def resource_path(relative_path):
    try:
        base_path = sys._MEIPASS
    except AttributeError:
        base_path = path.abspath(".")
    return path.join(base_path, relative_path)

logo_path = resource_path("assets/hermes_logo.png")
sys.setrecursionlimit(int(100000))

####################
# Literature Miner #
####################

def literature_miner(request_title, email, terms, keywords_temp, mandatory, nresults):

    # Processing User Input
    hermes_mode="(Default Mode)"
    if mandatory:
        hermes_mode="(Strict Mode)"
    request = request_title
    email = email
    terms = terms
    keywords = []
    keywords_temp = keywords_temp
    keywords_temp = keywords_temp.split(",")
    for i in keywords_temp:
        keywords.append(i.strip(" "))
    keywords_mandatory = bool(mandatory)
    retmax = 100000
    nresults = int(nresults)

    progress_label0["text"] = "Initiating Literature Miner v5..."
    window.update()

    # Identification to Entrez Server
    Entrez.email = email

    # Creating Timestamp
    now = datetime.now()
    curr_year = int(now.strftime("%Y"))
    date = datetime.today().date()

    # Creating Result Folder
    working_directory = my_dir.get()
    progress_label0["text"] = "Creating report folder..."
    window.update()
    af.create_reports_folder(my_dir.get(), request, date)
    report_dir = path.join("{working_directory}".format(working_directory = working_directory), "{date}_{request_id}".format(request_id = request, date = date))
    figures_dir = path.join("{report_dir}".format(report_dir = report_dir), "figures")
    results_dir = path.join("{report_dir}".format(report_dir = report_dir), "results")
    log_file = "{report_dir}/{request}_log.txt".format(report_dir = report_dir, request = request)

    # Creating Report Log
    with open(log_file, 'w') as f:
        f.write("Request Title: {request_title}\n".format(request_title = request))
        f.write("Query Terms: {terms}\n".format(terms = terms))
        f.write("Parsing Keywords: {keywords}\n".format(keywords = keywords))
        f.write("Mandatory Keywords: {mandatory}\n".format(mandatory = keywords_mandatory))
        f.write("Retmax: {retmax}\n".format(retmax = retmax))
        f.write("N Results: {nresults}\n".format(nresults = nresults))
        f.write("\nLog:\n")

    # Querying Pubmed PMC Database for PMCIDs
    progress_label0["text"] = "Querying Pubmed PMC database..."
    with open(log_file, 'a') as f:
        f.write("Querying Pubmed PMC database...\n")
    window.update()
    handle = Entrez.esearch(db = "pmc", term = terms, retmax = str(retmax), sort = "relevance")
    query_result = Entrez.read(handle)
    id_list = query_result["IdList"]
    hits = query_result["Count"]
    progress_label0["text"] = "Retrieving complete PMC entries for {hits} hits...".format(hits = int(hits))
    with open(log_file, 'a') as f:
        f.write("Retrieving complete PMC entries for {hits} hits...\n".format(hits = int(hits)))
        f.write("Hits: {id_list}\n".format(id_list = id_list))
    window.update()

    # Implementing Rate-limited Fetch (3 requests/second Max)
    def fetch_record(pmcid):
        try:
            time.sleep(random.uniform(0.3, 0.5))
            handle = Entrez.efetch(db="pmc", id=pmcid, retmode = "xml", rettype = "full")
            record = str(list(Medline.parse(handle)))
            return record
        except Exception as e:
            return None
    
    # Initiating Result Dataframe
    stats_df = pd.DataFrame(columns = ["PMCID", "PMID", "Title", "Year", "Associated Keywords", "Citations", "Keywords Count", "Score"])

    # Mining Full PMC Article
    def process_pmcid(pmcid):

        # Fetch Record
        try:
            record = fetch_record(pmcid)
            if record is None:
                with open(log_file, 'a') as f:
                    f.write("Entrez fetch failed for {pmcid}\n".format(pmcid = pmcid))
                return {"error": f"{pmcid}"}
        except:
            pass
        
        # Parsing Record with Beautiful Soup
        try:    
            bs_record = bs(record, "lxml")
            article_title = bs_record.find("article-title").get_text()
            journal_title = bs_record.find("journal-title").get_text()
            associated_keywords = [i.get_text().lower() for i in bs_record.find_all("kwd")]
            pub_year = bs_record.find("pub-date").find("year").get_text()
        except:
            with open(log_file, 'a') as f:
                f.write("Issue with PMCID {pmcid} data request: Could not parse record\n".format(pmcid = pmcid))
            pass

        # Getting Abstract
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

        # Getting PMID
        article_ids = bs_record.find_all('article-id')
        pmid = None
        doi = None
        for id in article_ids:
            if "pmid" in str(id):
                pmid = id.get_text()
            if "doi" in str(id):
                doi = id.get_text()

        # Getting Authors
        contribs = bs_record.find_all("contrib")
        authors = []
        try:
            for author in contribs:
                authors.append("{surname}, {name}".format(surname = author.find("name").find("surname").get_text(), name = author.find("name").find("given-names").get_text()))
        except:
            authors = ["NA"]

        # Getting Citations
        n_citations = af.get_citations(pmid)

        # Getting Main Text
        introduction_pos = str(bs_record).find("<title>Introduction</title>")
        if introduction_pos == -1:
            introduction_pos = str(bs_record).lower().find("<title>Introduction".lower())
        references_pos = str(bs_record).find("<title>References</title>")
        if references_pos == -1:
            references_pos = str(bs_record).lower().find("<title>References".lower())
        article = bs(str(bs_record)[introduction_pos:references_pos], "lxml").get_text().lower()

        # Counting Keywords Occurence in Main Text
        keywords_count = af.keyword_counter(keywords, article)

        # Scoring Article
        if n_citations is None or n_citations == "NA":
            n_citations = 0
        if keywords_mandatory:
            keywords_score = 1
            for count in keywords_count.values():
                keywords_score = keywords_score * (count/len(keywords_count.keys()))
            score = math.log2((((n_citations + 1) * (keywords_score))/(((curr_year - int(pub_year))*0.25 + 1))*10) + 1)
        else:
            score = math.log2((((n_citations + 1) * (sum(keywords_count.values())/len(keywords_count.keys()) + 1))/(((curr_year - int(pub_year))*0.25 + 1))*10) + 1)
        with open(log_file, 'a') as f:
            f.write("Processing complete for {pmcid}\n".format(pmcid = pmcid))

        # Return Article Mining Results
        return {"PMCID": pmcid,
                "PMID": pmid,
                "Title": article_title,
                "Year": int(pub_year),
                "Journal": journal_title,
                "DOI": doi,
                "Authors": authors,
                "Associated Keywords": associated_keywords,
                "Abstract": abstract,
                "Citations": n_citations,
                "Keywords Count": keywords_count,
                "Score": score}
    
    # Initiating Mining Result List
    results = []

    # Initiating Data Request Errors List
    errors = []
    
    # Multi-treaded Executor
    n_parsed = 0
    progress_bar.configure(mode='determinate', maximum=len(id_list), value=0)
    progress_label["text"] = "Mining in progress... [0/{total}]".format(total = len(id_list))
    progress_bar.update_idletasks()
    progress_label.update_idletasks()

    with ThreadPoolExecutor(10) as executor:
        futures = {executor.submit(process_pmcid, pmcid): pmcid for pmcid in id_list}
        for future in as_completed(futures):
            result = future.result()
            if "error" in result:
                errors.append(result["error"])
            else:
                results.append(result)
                n_parsed += 1
                progress_bar["value"] = n_parsed
                progress_label["text"] = "Mining in progress... [{n_parsed}/{total}]".format(n_parsed = n_parsed, total = len(id_list))
                progress_bar.update_idletasks()
                progress_label.update_idletasks()
    
    # Initiating Retry Counter
    retry_count = 0

    # Reattempt Miner if Data Fetch Error
    while (len(errors) > 0) & (retry_count < 10):

        # Initiating Temporary Error List
        errors_temp = []

        # Multi-treaded Executor
        with ThreadPoolExecutor(5) as executor:
            futures = {executor.submit(process_pmcid, pmcid): pmcid for pmcid in errors}
            for future in as_completed(futures):
                result = future.result()
                if "error" in result:
                    errors_temp.append(result["error"])
                else:
                    results.append(result)
                    n_parsed += 1
                    progress_bar["value"] = n_parsed
                    progress_label["text"] = "Mining in progress... [{n_parsed}/{total}]".format(n_parsed = n_parsed, total = len(id_list))
                    progress_bar.update_idletasks()
                    progress_label.update_idletasks()   
        errors = errors_temp
        retry_count += 1

    # Generating Result Dataframe
    results_df = pd.DataFrame(results)

    progress_bar["value"] = 100
    window.update()
    with open(log_file, 'a') as f:
                f.write("Processing done\n")
                f.write("Processing errors: {errors}\n".format(errors = errors))

    # Generating statistics for query
    hits = int(hits)
    results_df = results_df.drop_duplicates(subset=["PMCID"])
    parsed = len(results_df)
    progress_label["text"] = "Mining completed! [{succes}/{hits}]".format(succes = parsed, hits = hits)
    with open(log_file, 'a') as f:
        f.write("Mining completed! [{succes}/{hits}]\n".format(succes = parsed, hits = hits))
    progress_bar["value"] = 100
    window.update()

    # Writing Results Dataframe to Disk
    combined_df = pd.concat([stats_df, results_df], ignore_index=True)
    combined_df.to_csv("{results_dir}/literature_results.csv".format(results_dir=results_dir), index=False)

    # Generating Graphs
    progress_label2["text"] = "Generating report graphs..."
    with open(log_file, 'a') as f:
        f.write("Generating report graphs...\n")
    window.update()
    af.graph_pub_years(combined_df, request, figures_dir)
    af.graph_associated_keywords(combined_df, request, figures_dir)
    af.graph_stats_summary(combined_df, request, figures_dir)

    # Keeping Top Results for Report
    results_df = results_df.sort_values(by = "Score", ascending=False).reset_index(drop = True)
    results_df = results_df[results_df["Score"] > 0]
    if len(results_df) > nresults:
        results_df = results_df.head(nresults)

    # Generating AI Summary for Top Results with Multi-treaded Executor
    summaries = []
    n_summarized = 0
    progress_bar.configure(mode='determinate', maximum=len(results_df), value=0)
    progress_label2["text"] = "Summarizing top results... [0/{total}]".format(total = len(results_df))
    progress_bar.update_idletasks()
    progress_label2.update_idletasks()

    with ThreadPoolExecutor(3) as executor:
        futures = {executor.submit(summarize.summarize_article, results_df.loc[i, "PMCID"]): results_df.loc[i, "PMCID"] for i in range(len(results_df))}
        for future in as_completed(futures):
            result = future.result()
            if "error" in result:
                summaries.append(result["error"])
            else:
                summaries.append(result)
                n_summarized += 1
                progress_bar["value"] = n_summarized
                progress_label2["text"] = "Summarizing top results... [{n_summarized}/{total}]".format(n_summarized = n_summarized, total = len(results_df))
                progress_bar.update_idletasks()
                progress_label2.update_idletasks()
    
    summaries_df = pd.DataFrame(summaries, columns=["PMCID", "Summary"])
    results_df = results_df.merge(summaries_df, on="PMCID", how="left")

    # Generating the PDF report
    class PDF_Report(FPDF):
        def __init__(self):
            super().__init__()
            self.font = "Arial"
            self.m = 10 
            self.pw = 210 - 2*self.m
            self.line_thickness = 0
            self.height = 7

        def header(self):
            self.set_font(self.font, '', 10)
            self.cell(w=(self.pw/3), h=self.height, txt="", border=self.line_thickness, ln=0, align='L')
            self.cell(w=(self.pw/3), h=self.height, txt='Hermes (v1.1) - Literature Mining Report', border=self.line_thickness, ln=0, align = 'C')
            self.cell(w=(self.pw/3), h=self.height, txt="{date}".format(date = date), border=self.line_thickness, ln=1, align = "R")
            self.image(logo_path, 10, 8, 35)
            self.ln(5)

        def footer(self):
            self.set_y(-15)
            self.set_font(self.font, '', 8)
            self.cell(w=(self.pw), h=self.height, txt='Page {page_number}'.format(page_number = self.page_no()), border=self.line_thickness, ln=1, align = 'C')

    ## Creating first page of the report (Title + Query Statistics)
    progress_label2["text"] = "Generating PDF report..."
    with open(log_file, 'a') as f:
        f.write("Generating PDF graphs...\n")
    window.update()
    pdf = PDF_Report()
    pdf.add_page()
    #pdf.ln(1)
    pdf.set_font(pdf.font, 'B', 24)
    pdf.cell(w=pdf.pw, h=25, txt="Literature Mining Report", border=pdf.line_thickness, ln=1, align='C')
    pdf.set_font(pdf.font, 'B', 12)
    pdf.cell(w=pdf.pw, h=8, txt="Query Statistics", border=pdf.line_thickness, ln=1, align='L')
    pdf.set_font(pdf.font, '', 10)
    pdf.cell(w=pdf.pw, h=5, txt="Query date: {date}".format(date = date), border=pdf.line_thickness, ln=1, align='L')
    pdf.multi_cell(w=pdf.pw, h=5, txt='Query terms: "{terms}"'.format(terms = terms), border=pdf.line_thickness, align='L')
    pdf.cell(w=pdf.pw, h=5, txt='Articles parsed/hits: {parsed}/{hits} ({percentage}%)'.format(parsed = int(parsed), hits = int(hits), percentage = "{:.2f}".format(parsed*100/hits)), border=pdf.line_thickness, ln=1, align='L')
    pdf.cell(w=pdf.pw, h=5, txt='Scoring keywords: {keywords} {mode}'.format(keywords = ", ".join(keywords), mode = hermes_mode), border=pdf.line_thickness, ln=1, align='L')
    pdf.ln(3)
    pdf.image("{figures_dir}/{request}_pubyears.png".format(figures_dir = figures_dir, request = request), w = pdf.pw)
    pdf.image("{figures_dir}/{request}_asskeywords.png".format(figures_dir = figures_dir, request = request), w = pdf.pw)
    pdf.add_page()
    pdf.ln(2)
    pdf.image("{figures_dir}/{request}_stats_summary.png".format(figures_dir = figures_dir, request = request), w = pdf.pw)

    ## Generating the top 25 result table
    pdf.ln(2)
    ### Generating the table header
    pdf.set_font(pdf.font, 'B', 10)
    pdf.cell(w=pdf.pw, h=5, txt="Query Results", border=pdf.line_thickness, ln=1, align='C')
    pdf.ln(1)
    pdf.set_font(pdf.font, 'B', 9)
    pdf.table_border = 1
    pdf.cell(w=(3*pdf.pw/30), h= 5, txt="PMID", border=pdf.table_border, ln=0, align='C')
    pdf.cell(w=(2*pdf.pw/30), h= 5, txt="Year", border=pdf.table_border, ln=0, align='C')
    pdf.cell(w=(4*pdf.pw/30), h= 5, txt="Journal", border=pdf.table_border, ln=0, align='C')
    pdf.cell(w=(5*pdf.pw/30), h= 5, txt="Authors", border=pdf.table_border, ln=0, align='C')
    pdf.cell(w=(11*pdf.pw/30), h= 5, txt="Title", border=pdf.table_border, ln=0, align='C')
    pdf.cell(w=(3*pdf.pw/30), h= 5, txt="Citations", border=pdf.table_border, ln=0, align='C')
    pdf.cell(w=(2*pdf.pw/30), h= 5, txt="Score", border=pdf.table_border, ln=1, align='C')
    ### Generating the table from the dataframe
    pdf.set_font(pdf.font, '', 8)
    for i in range(len(results_df)):
        pdf.cell(w=(3*pdf.pw/30), h= 5, txt=results_df["PMID"][i], link = "https://doi.org/{doi}".format(doi = results_df["DOI"][i]), border=pdf.table_border, ln=0, align='C')
        pdf.cell(w=(2*pdf.pw/30), h= 5, txt=str(results_df["Year"][i]), border=pdf.table_border, ln=0, align='C')
        pdf.cell(w=(4*pdf.pw/30), h= 5, txt=af.filter_string(results_df["Journal"][i][:15]), border=pdf.table_border, ln=0, align='L')
        if len(results_df["Authors"][i]) > 1:
            pdf.cell(w=(5*pdf.pw/30), h= 5, txt=af.filter_string("{first_author} et al.".format(first_author = results_df["Authors"][i][0][:15])), border=pdf.table_border, ln=0, align='L')
        else:
            if len(results_df["Authors"][i]) < 1:
                pdf.cell(w=(5*pdf.pw/30), h= 5, txt="NaN", border=pdf.table_border, ln=0, align='L')
            else:
                pdf.cell(w=(5*pdf.pw/30), h= 5, txt=af.filter_string(results_df["Authors"][i][0][:15]), border=pdf.table_border, ln=0, align='L')
        if len(results_df["Title"][i]) > 45:
            pdf.cell(w=(11*pdf.pw/30), h= 5, txt=af.filter_string("{title}...".format(title = results_df["Title"][i][:45].lower())), border=pdf.table_border, ln=0, align='L')
        else:
            pdf.cell(w=(11*pdf.pw/30), h= 5, txt=af.filter_string(results_df["Title"][i]), border=pdf.table_border, ln=0, align='L')
        pdf.cell(w=(3*pdf.pw/30), h= 5, txt="{citation}".format(citation = results_df["Citations"][i]), border=pdf.table_border, ln=0, align='C')
        pdf.cell(w=(2*pdf.pw/30), h= 5, txt="{score}".format(score = "{:.2f}".format(results_df["Score"][i])), border=pdf.table_border, ln=1, align='C')

    ## Generating long reports
    counter = 1
    pdf.add_page()
    for i in range(len(results_df)):
        if counter > 1:
            pdf.add_page()
            pdf.ln(3)
            counter = 1
        pdf.set_font(pdf.font, 'B', 10)
        pdf.cell(w=pdf.pw, h=5, txt="Query Result {i}/{len_results}".format(i = i + 1, len_results = len(results_df)), border=pdf.line_thickness, align='L')
        pdf.ln(5)
        pdf.set_font(pdf.font, '', 10)
        pdf.multi_cell(w=pdf.pw, h=5, txt=af.filter_string(results_df["Title"][i]), border=pdf.line_thickness, align='L')
        pdf.set_font(pdf.font, '', 9)
        pdf.multi_cell(w=pdf.pw, h=5, txt=af.filter_string("; ".join(results_df["Authors"][i])), border=pdf.line_thickness, align='L')
        pdf.set_font(pdf.font, '', 8)
        pdf.cell(w=pdf.pw, h=5, txt=af.filter_string(results_df["Journal"][i]), border=pdf.line_thickness, ln=1, align='L')
        pdf.cell(w=(4*pdf.pw/30), h= 5, txt="Year: {year}".format(year = results_df["Year"][i]), border=pdf.line_thickness, ln=0, align='L')
        pdf.cell(w=(4*pdf.pw/30), h= 5, txt="PMID: {pmid}".format(pmid = results_df["PMID"][i]), border=pdf.line_thickness, ln=0, align='L')
        pdf.cell(w=(8*pdf.pw/30), h= 5, txt="doi: {doi}".format(doi = results_df["DOI"][i]), link = "https://doi.org/{doi}".format(doi = results_df["DOI"][i]), border=pdf.line_thickness, ln=0, align='L')
        pdf.cell(w=(4*pdf.pw/30), h= 5, txt="Citations: {citations}".format(citations = results_df["Citations"][i]), border=pdf.line_thickness, ln=0, align='L')
        pdf.cell(w=(4*pdf.pw/30), h= 5, txt="Score: {score}".format(score = "{:.2f}".format(results_df["Score"][i])), border=pdf.line_thickness, ln=1, align='L')
        pdf.cell(w=(4*pdf.pw/30), h= 5, txt=af.filter_string("Keywords hits: {keywords}".format(keywords = results_df["Keywords Count"][i])), border=pdf.line_thickness, ln=1, align='L')
        pdf.ln(2)
        pdf.set_font(pdf.font, 'B', 8)
        pdf.cell(w=pdf.pw, h=5, txt="Abstract:", border=pdf.line_thickness, align='L')
        pdf.ln(4)
        pdf.set_font(pdf.font, '', 8)
        pdf.multi_cell(w=pdf.pw, h=5, txt=af.filter_string(results_df["Abstract"][i]), border=pdf.line_thickness, align='L')
        pdf.ln(2)
        pdf.set_font(pdf.font, 'B', 8)
        pdf.cell(w=pdf.pw, h=5, txt="Summary:", border=pdf.line_thickness, align='L')
        pdf.ln(4)
        pdf.set_font(pdf.font, '', 8)
        pdf.multi_cell(w=pdf.pw, h=5, txt=af.filter_string(results_df["Summary"][i][0]), border=pdf.line_thickness, align='L')
        pdf.ln(2)
        pdf.set_font(pdf.font, 'B', 8)
        pdf.cell(w=pdf.pw, h=5, txt="Figures:", border=pdf.line_thickness, align='L')
        pdf.ln(4)
        pdf.set_font(pdf.font, '', 8)
        for j in range(len(results_df["Summary"][i][9])):
            pdf.multi_cell(w=pdf.pw, h=5, txt="Figure {j}: ".format(j = j+1) + af.filter_string(results_df["Summary"][i][9][j]), border=pdf.line_thickness, align='L')
            pdf.ln(1)
        pdf.ln(1)
        pdf.set_font(pdf.font, 'B', 8)
        pdf.cell(w=pdf.pw, h=5, txt="Mentioned biomedical entities:", border=pdf.line_thickness, align='L')
        pdf.ln(4)
        pdf.set_font(pdf.font, '', 8)
        pdf.multi_cell(w=pdf.pw, h=5, txt="Genes: " + af.filter_string(results_df["Summary"][i][1]), border=pdf.line_thickness, align='L')
        pdf.ln(1)
        pdf.multi_cell(w=pdf.pw, h=5, txt="Proteins: " + af.filter_string(results_df["Summary"][i][2]), border=pdf.line_thickness, align='L')
        pdf.ln(1)
        pdf.multi_cell(w=pdf.pw, h=5, txt="Diseases: " + af.filter_string(results_df["Summary"][i][3]), border=pdf.line_thickness, align='L')
        pdf.ln(1)
        pdf.multi_cell(w=pdf.pw, h=5, txt="Chemicals: " + af.filter_string(results_df["Summary"][i][4]), border=pdf.line_thickness, align='L')
        pdf.ln(1)
        pdf.multi_cell(w=pdf.pw, h=5, txt="Cells: " + af.filter_string(results_df["Summary"][i][5]), border=pdf.line_thickness, align='L')
        pdf.ln(1)
        pdf.multi_cell(w=pdf.pw, h=5, txt="Organisms: " + af.filter_string(results_df["Summary"][i][6]), border=pdf.line_thickness, align='L')
        pdf.ln(1)
        pdf.multi_cell(w=pdf.pw, h=5, txt="Tissues: " + af.filter_string(results_df["Summary"][i][7]), border=pdf.line_thickness, align='L')
        pdf.ln(1)
        pdf.multi_cell(w=pdf.pw, h=5, txt="Pathways: " + af.filter_string(results_df["Summary"][i][8]), border=pdf.line_thickness, align='L')
        counter += 1

    # Creating the output PDF report
    pdf.output("{report_dir}/{date}_{request}.pdf".format(report_dir = report_dir, date = date, request = request), 'F')
    submit_label["text"] = ""
    progress_label2["text"] = "Report ready!"
    with open(log_file, 'a') as f:
        f.write("Report ready!\n")
        f.write(":)")
    progress_bar["value"] = 0
    window.update() 

###########################################################################################################################################
# GUI definition
window = tk.Tk()
window.minsize(1000, 650)
window.configure(bg='white')
window.title("Hermes (v1.1) - Open Source Literature Mining")
style = ttk.Style(window)
style.theme_use('default')
style.configure("custom.Horizontal.TProgressbar", troughcolor='white', background='#1F4063', thickness=20)
default_font = tkFont.nametofont("TkDefaultFont")
default_font.configure(family="Arial", size=9)

# Frame
frame = tk.Frame(window, bg="white")
frame.pack(padx=10, pady=10, anchor="nw")

# Header
img_logo = tk.PhotoImage(file=logo_path)
logo = tk.Label(master = frame, image=img_logo, bg='white', fg='white')
logo.grid(row=0, column=0, sticky="w")
greeting = tk.Label(master = frame, text = "Hermes is an open-source mining tool for open-access literature, enabling scoring, and ranking of full-text articles based on customizable keyword sets and other relevance metrics.", bg='white', fg='black')
greeting.grid(row=2, column=0, sticky="w")

# Request Title Fields
request_title = tk.StringVar(master = window, value = None)
request_title_ok = tk.BooleanVar(master = window, value = False)
title_text = tk.Label(master = frame, text = "Enter request title:", bg='white', fg='black')
title_text.grid(row=3, column=0, sticky="w", pady=(15,0))
greeting_title_1 = tk.Label(master = frame, text = "Specify a unique request title for report generation (e.g. 't.reesei_cellulases'). Reusing a title will overwrite the existing report.", bg='white', fg='black')
greeting_title_1.grid(row=4, column=0, sticky="w")
title_frame = tk.Frame(frame, bg="white")
title_frame.grid(row=5, column=0, sticky="w")
title_entry = tk.Entry(master = title_frame, relief='flat', width = 50, highlightthickness=1, bd=0, highlightbackground="black")
title_entry.grid(row=1, column=0, sticky="w", padx=3)
title_entry_button = tk.Button(master = title_frame, text = ">>>", command = title_request_button, width=3, height=1, bg='white', fg="black", relief='raised')
title_entry_button.grid(row=1, column=1, sticky="w", padx=5)
title_entry_field = tk.Label(master = title_frame, text = "[Request Title]", bg='white', fg='black')
title_entry_field.grid(row=1, column=2, sticky="w", padx=2)

## Identification to Entrez server. PubMed requires an email address (string) associated with each query.
email = tk.StringVar(master = window, value = None)
email_ok = tk.BooleanVar(master = window, value = False)
email_text = tk.Label(master = frame, text = "Enter email address:", bg='white', fg='black')
email_text.grid(row=6, column=0, sticky="w", pady=(10,0))
email_text_2 = tk.Label(master = frame, text = "Enter your email address for identification with the Entrez server.", bg='white', fg='black')
email_text_2.grid(row=7, column=0, sticky="w")
email_frame = tk.Frame(frame, bg="white")
email_frame.grid(row=8, column=0, sticky="w")
email_entry = tk.Entry(master = email_frame, relief='flat', width = 50, highlightthickness=1, bd=0, highlightbackground="black")
email_entry.grid(row=1, column=0, sticky="w", padx=3)
email_entry_button = tk.Button(master = email_frame, text = ">>>", command = email_button, width=3, height=1, bg='white', fg="black", relief='raised')
email_entry_button.grid(row=1, column=1, sticky="w", padx=5)
email_entry_field = tk.Label(master = email_frame, text = "[Email]", bg='white', fg='black')
email_entry_field.grid(row=1, column=2, sticky="w", padx=2)

## Defining Pubmed research terms (string). Can include PubMed filters to refine the search (https://pubmed.ncbi.nlm.nih.gov/help/#help-filters).
terms = tk.StringVar(master = window, value = None)
terms_ok = tk.BooleanVar(master = window, value = False)
terms_text = tk.Label(master = frame, text = "Enter query terms:", bg='white', fg='black')
terms_text.grid(row=9, column=0, sticky="w", pady=(10,0))
terms_text_2 = tk.Label(master = frame, text = "Define your Pubmed query (e.g. 'cellulase production trichoderma reesei'). PubMed filters are supported (https://pubmed.ncbi.nlm.nih.gov/help/#help-filters).", bg='white', fg='black')
terms_text_2.grid(row=10, column=0, sticky="w")
terms_frame = tk.Frame(frame, bg="white")
terms_frame.grid(row=11, column=0, sticky="w")
terms_entry = tk.Entry(master = terms_frame, width = 50, relief='flat', highlightthickness=1, bd=0, highlightbackground="black")
terms_entry.grid(row=1, column=0, sticky="w", padx=3)
terms_entry_button = tk.Button(master = terms_frame, text = ">>>", command = terms_button, width=3, height=1, bg='white', fg="black", relief='raised')
terms_entry_button.grid(row=1, column=1, sticky="w", padx=5)
terms_entry_field = tk.Label(master = terms_frame, text = "[Query Terms]", bg='white', fg='black')
terms_entry_field.grid(row=1, column=2, sticky="w", padx=2)

## Defining keywords (list of strings) for parsing and scoring. Keywords mentadory (True/False).
keywords_temp = tk.StringVar(master = window, value = None)
keywords_ok = tk.BooleanVar(master = window, value = False)
keywords_text = tk.Label(master = frame, text = "Enter scoring keywords:", bg='white', fg='black')
keywords_text.grid(row=12, column=0, sticky="w", pady=(10,0))
keywords_text_2 = tk.Label(master = frame, text = "Define keywords for article scoring (e.g. 'cellulase, reesei, cbh1, xyr1').", bg='white', fg='black')
keywords_text_2.grid(row=13, column=0, sticky="w")
keywords_frame = tk.Frame(frame, bg="white")
keywords_frame.grid(row=14, column=0, sticky="w")
keywords_entry = tk.Entry(master = keywords_frame, width = 50, relief='flat', highlightthickness=1, bd=0, highlightbackground="black")
keywords_entry.grid(row=1, column=0, sticky="w", padx=3)
keywords_entry_button = tk.Button(master = keywords_frame, text = ">>>", command = keywords_button, width=3, height=1, bg='white', fg="black", relief='raised')
keywords_entry_button.grid(row=1, column=1, sticky="w", padx=5)
keywords_entry_field = tk.Label(master = keywords_frame, text = "[Scoring Keywords]", bg='white', fg='black')
keywords_entry_field.grid(row=1, column=2, sticky="w", padx=2)
keywords_mandatory = tk.BooleanVar(master = window, value = False)
keywords_check = tk.IntVar()
keywords_check_button = tk.Checkbutton(master = frame, text = "Strict Mode (All keywords required in article main text)", variable = keywords_check, onvalue = 1, offvalue = 0, command = mandatory_keywords_check, bg='white', fg="black")
keywords_check_button.grid(row=15, column=0, sticky="w")

## Defining the number of results to include in the report.
nresults = tk.IntVar(master = window, value = 25)
nresults_ok = tk.BooleanVar(master = window, value = True)
nresults_text_0 = tk.Label(master = frame, text = "Enter number of results:", bg='white', fg='black')
nresults_text_0.grid(row=16, column=0, sticky="w", pady=(10,0))
nresults_text = tk.Label(master = frame, text = "Define a number of articles to be included in the report.", bg='white', fg='black')
nresults_text.grid(row=17, column=0, sticky="w")
nresults_frame = tk.Frame(frame, bg="white")
nresults_frame.grid(row=18, column=0, sticky="w")
nresults_entry = tk.Entry(master = nresults_frame, width = 5, relief='flat', highlightthickness=1, bd=0, highlightbackground="black")
nresults_entry.grid(row=1, column=0, sticky="w", padx=3)
nresults_entry_button = tk.Button(master = nresults_frame, text = ">>>", command = nresults_button, width=3, height=1, bg='white', fg="black", relief='raised')
nresults_entry_button.grid(row=1, column=1, sticky="w", padx=5)
nresults_entry_field = tk.Label(master = nresults_frame, text = str(nresults.get()), bg='white', fg='black')
nresults_entry_field.grid(row=1, column=2, sticky="w", padx=2)

# Submission
submit_frame = tk.Frame(frame, bg="white")
submit_frame.grid(row=19, column=0, sticky="w", pady=(15,0))
submit_button = tk.Button(master = submit_frame, text = "Submit Query", command = submit_query, bg='white', fg="black", relief='raised')
submit_button.grid(row=1, column=0, sticky="w", padx=2)
submit_label = tk.Label(master = submit_frame, text = "", bg='white', fg="black")
submit_label.grid(row=1, column=1, sticky="w", padx=3)


# Report status
progress_bar = ttk.Progressbar(master = frame, orient = tk.HORIZONTAL, style="custom.Horizontal.TProgressbar", mode = "determinate", length = 1000)
progress_bar.grid(row=20, column=0, sticky="w", padx=2, pady=10)
progress_label0 = tk.Label(master = frame, text = "", bg='white', fg='black')
progress_label0.grid(row=21, column=0, sticky="w")
progress_label = tk.Label(master = frame, text = "", bg='white', fg='black')
progress_label.grid(row=22, column=0, sticky="w")
progress_label2 = tk.Label(master = frame, text = "", bg='white', fg='black')
progress_label2.grid(row=23, column=0, sticky="w")

# Save directory
my_dir = tk.StringVar(master = window, value = "")

window.mainloop()