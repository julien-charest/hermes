#!/usr/bin/env python3

####################################################################
# Hermes v1.2 - Open-source mining tool for open-access literature #
# 2025-12-18                                                       #
# Written by Julien Charest & Katarina Priselac                    #
####################################################################

from os import path, getcwd, getpid
import sys
import ssl
ssl._create_default_https_context = ssl._create_unverified_context
import tkinter as tk
from tkinter import ttk, filedialog
import tkinter.font as tkFont
from datetime import datetime
from Bio import Entrez
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed
import time
import psutil


#########################
# Import HERMES Modules #
#########################
import modules.adding_functions as af
from modules.mining import mine_pmcid
from modules.generating_report import generate_pdf_report
import modules.summarizing as summarizing
from modules.generating_graphs import *
from modules.scoring import *


#########################
# Configure NCBI Access #
#########################

# Without an API key (e.g. "64b7027d492e3d65e5d614cc027c7c247408"), Entrez requests are limited to ~3 per second. (default = None)
# Supplying a valid API key increases the rate limit to ~10 requests/second.
Entrez.api_key = None 

# Maximum number of records returned per Entrez query (default = 100000)
retmax = 100000


##################################
# Configure Performance Settings #
##################################

# Parallel execution settings
max_workers_mining = 10     # Number of concurrent workers for mining tasks.
max_workers_summary = 10    # Number of concurrent workers for summarization tasks.

# Summarization configuration in modules/summarizing.py


################################
# Global Runtime Configuration #
################################

# Resolve path to bundled logo
logo_path=af.ressource_path("assets/hermes_logo.png")

# Increase recursion limit to accommodate deep parsing and processing steps
sys.setrecursionlimit(int(100000))

# Initialize process handle and peak resident memory tracker
proc = psutil.Process(getpid())
peak_rss = 0


#####################################
# Taking User Input from HERMES App #
#####################################

def title_request_button():
    
    """Validate the request title entered by the user and update the GUI state."""

    request_title.set(str(title_entry.get()))
    if len(request_title.get()) > 0:
        request_title_ok.set(True)
        title_entry_field["text"] = request_title.get()
    else:
        request_title_ok.set(False)
        title_entry_field["text"] = "Invalid input. Please define a request title."


def email_button():

    """Validate the user-provided email address required for Entrez requests and update the GUI state."""

    email.set(str(email_entry.get()))
    if all(x in email.get() for x in ["@", "."]):
        email_ok.set(True)
        email_entry_field["text"] = email.get()
    else:
        email_ok.set(False)
        email_entry_field["text"] = "Invalid input. Please enter a valid email address: "


def terms_button():

    """Validate PubMed Central query terms entered by the user and update the GUI state."""

    terms.set(str(terms_entry.get()))
    if type(terms.get()) == str and len(terms.get()) > 0:
        terms_ok.set(True)
        terms_entry_field["text"] = terms.get()
    else:
        terms_ok.set(False)
        terms_entry_field["text"] = "Invalid input. Please define PubMed PMC query terms."


def keywords_button():

    """Parse and validate a comma-separated list of keywords entered by the user and update the GUI state."""

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
                keywords_list_text += "{keyword}, ".format(keyword=keywords[j])
            keywords_list_text += keywords[-1]
            keywords_temp.set(keywords_list_text)
            keywords_entry_field["text"] = keywords_temp.get()
    else:
        keywords_temp.set("")
        keywords_ok.set(False)
        keywords_entry_field["text"] = "Invalid input. Please define a valid keywords list."


def mandatory_keywords_check():

    """Toggle strict keyword matching mode based on the user's selection."""

    if keywords_check.get() == 1:
        keywords_mandatory.set(True)
    else:
        keywords_mandatory.set(False)


def nresults_button():

    """Validate the number of results to include in the final report and update the GUI state."""

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


def submit_query():

    """Validate all user inputs and, if valid, launch the HERMES pipeline and close the GUI."""

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


def update_peak_memory():

    """Update the tracked peak resident memory (RSS) for the current process."""

    global peak_rss
    rss = proc.memory_info().rss
    peak_rss = max(peak_rss, rss)


def launch_script():

    """Collect validated user inputs from the GUI and initiate the HERMES literature mining pipeline."""

    request_var = request_title.get()
    email_var = email.get()
    terms_var = terms.get()
    keywords_temp_var = keywords_temp.get()
    mandatory_var = keywords_mandatory.get()
    nresults_var = nresults.get()
    print("[INFO] Initializing HERMES v1.2...")
    literature_miner(request_var, email_var, terms_var, keywords_temp_var, mandatory_var, nresults_var, retmax)


####################################
# HERMES — Main Execution Pipeline #
####################################

def literature_miner(request, email, terms, keywords_temp, mandatory, nresults, retmax):

    """Execute the HERMES literature mining workflow from query submission to report generation."""

    # Record start time for runtime measurement
    t0 = time.perf_counter()

    # Update peak memory usage
    update_peak_memory()

    # Normalize and validate user-provided inputs and execution mode
    print("[STATUS] Processing user input...")

    # Determine HERMES execution mode
    hermes_mode = "(Strict Mode)" if mandatory else "(Default Mode)"
    if mandatory:
        print("[INFO] Strict mode enabled.")

    # Normalize inputs / parameters
    keywords = [k.strip() for k in str(keywords_temp).split(",") if k.strip()]
    keywords_mandatory = bool(mandatory)
    nresults = int(nresults)
    retmax = int(retmax)

    # Update GUI
    progress_label0["text"] = "Initializing HERMES v1.2..."
    window.update()

    # Configure NCBI Entrez access (user identification and rate limiting)
    print("[STATUS] Connecting to Entrez server...")

    Entrez.email = email
    if Entrez.api_key:
        print("[INFO] Using NCBI API key: rate limit ~10 requests/second.")
    else:
        print("[INFO] No NCBI API key detected: rate limit ~3 requests/second.")

    # Record current date for report metadata and output labeling
    date = datetime.today().date()

    # Update GUI
    print("[STATUS] Creating report directory...")
    progress_label0["text"] = "Creating report directory..."
    window.update()

    # Initialize report output directories (base, figures, and results)
    working_directory = my_dir.get()
    af.create_reports_folder(working_directory, request, date)
    report_dir = path.join(working_directory, f"{date}_{request}")
    figures_dir = path.join(report_dir, "figures")
    results_dir = path.join(report_dir, "results")

    # Initialize run-specific log file capturing query parameters and execution context
    print("[STATUS] Creating log file...")
    log_file = path.join(report_dir, f"{request}_log.txt")
    with open(log_file, "w", encoding="utf-8") as f:
        f.write(f"Request Title: {request}\n")
        f.write(f"Query Terms: {terms}\n")
        f.write(f"Scoring Keywords: {keywords}\n")
        f.write(f"Strict Mode Enabled: {keywords_mandatory}\n")
        f.write(f"Retmax: {retmax}\n")
        f.write(f"N Results: {nresults}\n")
        f.write("\nLog:\n")

    # Update peak memory usage
    update_peak_memory()

    # Update GUI
    print("[STATUS] Querying PubMed Central database...")
    progress_label0["text"] = "Querying PubMed Central database..."
    window.update()

    # Update log
    with open(log_file, "a", encoding="utf-8") as f:
        f.write("[STATUS] Querying PubMed Central database...\n")

    # Query PubMed Central via Entrez to retrieve matching PMC identifiers
    handle = Entrez.esearch(db = "pmc", term = terms, retmax = str(retmax), sort = "relevance")
    query_result = Entrez.read(handle)
    id_list = query_result["IdList"]
    hits = query_result["Count"]

    # Update GUI
    print("[INFO] Retrieving complete PMC entries for {hits} hits...".format(hits = int(hits)))
    progress_label0["text"] = "Retrieving complete PMC entries for {hits} hits...".format(hits = int(hits))
    window.update()
    
    # Update log
    with open(log_file, 'a') as f:
        f.write("[INFO] Retrieving complete PMC entries for {hits} hits...\n".format(hits = int(hits)))
        f.write("[INFO] Hits: {id_list}\n".format(id_list = id_list))
    
    # Initialize containers for mining results and request-level errors
    results = []
    errors = []
    
    # Parallelize PMC article mining with a thread pool while tracking progress and logging per-article failures
    print("[STATUS] Literature mining in progress...")
    n_parsed = 0
    progress_bar.configure(mode='determinate', maximum=len(id_list), value=0)
    progress_label["text"] = "Literature mining in progress... [0/{total}]".format(total = len(id_list))
    progress_bar.update_idletasks()
    progress_label.update_idletasks()

    # Record start times for CPU time and wall-clock time measurement
    t1_pt = time.process_time()
    t1_pc = time.perf_counter()
    
    # Update peak memory usage
    update_peak_memory()

    ## Submit all PMC mining tasks to the executor and process results asynchronously as they complete
    with ThreadPoolExecutor(max_workers_mining) as executor:
        futures = {executor.submit(mine_pmcid, pmcid, keywords): pmcid for pmcid in id_list}
        for future in as_completed(futures):
            pmcid = futures[future]
            try:
                result = future.result()
            except Exception as e:
                errors.append(pmcid)
                with open(log_file, 'a') as f:
                    f.write(f"[ERROR] while processing {pmcid}: {e}\n")
                continue

            results.append(result)
            n_parsed += 1

            # Update GUI
            progress_bar["value"] = n_parsed
            progress_label["text"] = "Literature mining in progress... [{n_parsed}/{total}]".format(n_parsed=n_parsed, total=len(id_list))
            progress_bar.update_idletasks()
            progress_label.update_idletasks()
            
            # Update peak memory usage
            update_peak_memory()

    # Retry failed PMC mining tasks with a bounded number of attempts to mitigate transient fetch/parsing errors
    retry_count = 0

    ## Reattempt mining of failed PMC entries until all succeed or retry limit is reached
    while (len(errors) > 0) & (retry_count < 10):

        # Track PMC identifiers that fail again during the current retry iteration
        errors_temp = []

        # Re-dispatch failed PMC mining tasks in parallel (reduced to 5 workers)
        with ThreadPoolExecutor(5) as executor:
            futures = {executor.submit(mine_pmcid, pmcid, keywords): pmcid for pmcid in errors}
            for future in as_completed(futures):
                pmcid = futures[future]
                try:
                    result = future.result()
                except Exception as e:
                    errors_temp.append(pmcid)
                    with open(log_file, 'a') as f:
                        f.write(f"[ERROR] while processing {pmcid}: {e}\n")
                    continue

                results.append(result)
                n_parsed += 1

                # Update GUI
                progress_bar["value"] = n_parsed
                progress_label["text"] = "Literature mining in progress... [{n_parsed}/{total}]".format(n_parsed=n_parsed, total=len(id_list))
                progress_bar.update_idletasks()
                progress_label.update_idletasks()
                
                # Update peak memory usage
                update_peak_memory()
 
        errors = errors_temp
        retry_count += 1

    # Finalize mining stage and report cumulative processing time and coverage
    elapsed_pc = time.perf_counter() - t1_pc
    elapsed_pt = time.process_time() - t1_pt
    print("[STATUS] Literature mining complete.")
    print("[INFO] Processed {n_parsed}/{total} entries in {elapsed_pc:.2f} s (CPU time: {elapsed_pt:.2f} s).".format(n_parsed=n_parsed, total=len(id_list), elapsed_pc=elapsed_pc, elapsed_pt=elapsed_pt))

    # Update peak memory usage
    update_peak_memory()
    print(f"[INFO] Peak memory usage: {peak_rss / 1024**2:.1f} MiB")

    # Update log
    with open(log_file, 'a') as f:
        f.write("[STATUS] Literature mining complete.\n")
        f.write("[INFO] Processed {n_parsed}/{total} entries in {elapsed_pc:.2f} s (CPU time: {elapsed_pt:.2f} s).\n".format(n_parsed=n_parsed, total=len(id_list), elapsed_pc=elapsed_pc, elapsed_pt=elapsed_pt))
        f.write("[INFO] Peak memory usage: {mem:.1f} MiB\n".format(mem=peak_rss / 1024**2))

    # Consolidate mined article metadata into a structured results table for downstream analysis
    results_df = pd.DataFrame(results)

    # Update GUI
    progress_bar["value"] = 100
    window.update()

    # Record any PMCIDs that could not be processed after retries
    if len(errors) > 0:
        with open(log_file, 'a') as f:
            f.write("[ERROR] Articles not processed: {errors}\n".format(errors = errors))

    # Compute summary statistics and remove duplicate entries prior to scoring
    hits = int(hits)
    results_df = results_df.drop_duplicates(subset=["PMCID"])
    parsed = len(results_df)

    # Update GUI
    progress_label["text"] = "Literature mining complete. [{succes}/{hits}]".format(succes = parsed, hits = hits)
    progress_bar["value"] = 100
    window.update()

    # Score and rank retrieved articles based on keyword relevance (IDF-weighted)
    results_df, idf = scoring_articles(results_df, keywords_mandatory)
    
    # Update peak memory usage
    update_peak_memory()

    # Persist finalized mining results to disk for downstream inspection and reuse
    print("[STATUS] Writing mining results to disk...")

    with open(log_file, 'a') as f:
        f.write("[STATUS] Writing mining results to disk...\n")

    results_df.to_csv("{results_dir}/hermes_results.csv".format(results_dir=results_dir), index=False)

    # Generate summary visualizations for the analysis report
    print("[STATUS] Generating report graphs...")
    progress_label2["text"] = "Generating report graphs..."
    window.update()

    with open(log_file, 'a') as f:
        f.write("[STATUS] Generating report graphs...\n")

    graph_pub_years(results_df, request, figures_dir)
    graph_associated_keywords(results_df, request, figures_dir)
    graph_stats_summary(results_df, request, figures_dir)

    # Select and curate the top-scoring articles for inclusion in the final report
    results_df = results_df.sort_values(by = "Score", ascending=False).reset_index(drop = True)
    results_df = results_df[results_df["Score"] > 0]
    if len(results_df) > nresults:
        results_df = results_df.head(nresults)

    # Update peak memory usage
    update_peak_memory()

    # Record CPU and wall-clock start times for performance monitoring
    t2_pc = time.perf_counter()
    t2_pt = time.process_time()

    # Generate AI-based summaries for the top-ranked articles using concurrent execution
    errors = []
    top_results = results_df["PMCID"].tolist()
    summaries = {}
    n_summarized = 0

    print("[STATUS] Preparing summaries of {nresults} top results...".format(nresults = len(top_results)))

    with open(log_file, 'a') as f:
        f.write("[STATUS] Preparing summaries of {nresults} top results...\n".format(nresults = len(top_results)))

    progress_bar.configure(mode='determinate', maximum=len(top_results), value=0)
    progress_label2["text"] = "Preparing summaries of top results... [0/{total}]".format(total = len(top_results))
    progress_bar.update_idletasks()
    progress_label2.update_idletasks()

    ## Dispatch summarization tasks concurrently to improve throughput for independent articles
    with ThreadPoolExecutor(max_workers_summary) as executor:
        futures = {executor.submit(summarizing.summarize_article, pmcid): pmcid for pmcid in top_results}

        for future in as_completed(futures):
            pmcid = futures[future]
            try:
                result = future.result()
                summaries[pmcid] = result
            except Exception as e:
                errors.append(pmcid)
                with open(log_file, 'a') as f:
                    f.write(f"[ERROR] while summarizing {pmcid}: {e}\n")
                continue

            results.append(result)
            n_summarized += 1

            # Update GUI
            progress_bar["value"] = n_summarized
            progress_label2["text"] = "Preparing summaries of top results... [{n_summarized}/{total}]".format(n_summarized = n_summarized, total = len(results_df))
            progress_bar.update_idletasks()
            progress_label2.update_idletasks()

            # Update peak memory usage
            update_peak_memory()

    ## Retry failed summarizations with a bounded number of attempts to mitigate transient model/runtime errors
    retry_count = 0
    max_retries = 3

    while (len(errors) > 0) and (retry_count < max_retries):

        errors_temp = []

        # Re-dispatch failed summarization tasks in parallel
        with ThreadPoolExecutor(5) as executor:
            futures = {executor.submit(summarizing.summarize_article, pmcid): pmcid for pmcid in errors}

            for future in as_completed(futures):
                pmcid = futures[future]
                try:
                    result = future.result()
                    summaries[pmcid] = result
                except Exception as e:
                    errors_temp.append(pmcid)
                    with open(log_file, 'a') as f:
                        f.write(f"[ERROR] while summarizing {pmcid} (retry): {e}\n")
                    continue

                n_summarized += 1

                # Update GUI
                progress_bar["value"] = n_summarized
                progress_label2["text"] = "Preparing summaries of top results... [{n_summarized}/{total}]".format(n_summarized=n_summarized, total=len(top_results))
                progress_bar.update_idletasks()
                progress_label2.update_idletasks()

                # Update peak memory usage
                update_peak_memory()

        errors = errors_temp
        retry_count += 1

    # Finalize summarization stage: record runtime and integrate summaries into the results table
    elapsed_pc = time.perf_counter() - t2_pc
    elapsed_pt = time.process_time() - t2_pt
    print("[STATUS] Summarization of top-ranked articles complete.")
    print("[INFO] Summarized {n_summarized}/{total} entries in {elapsed_pc:.2f} s (CPU time: {elapsed_pt:.2f} s).".format(n_summarized = n_summarized, total = len(results_df), elapsed_pc=elapsed_pc, elapsed_pt=elapsed_pt))

    # Update peak memory usage
    update_peak_memory()
    print(f"[INFO] Peak memory usage: {peak_rss / 1024**2:.1f} MiB")

    # Update log
    with open(log_file, 'a') as f:
                f.write("[STATUS] AI summarization of top-ranked articles complete.\n")
                f.write("[INFO] Summarized {n_summarized}/{total} entries in {elapsed_pc:.2f} s (CPU time: {elapsed_pt:.2f} s).\n".format(n_summarized = n_summarized, total = len(results_df), elapsed_pc=elapsed_pc, elapsed_pt=elapsed_pt))
                f.write(f"[INFO] Peak memory usage: {peak_rss / 1024**2:.1f} MiB\n")

    # Update GUI
    results_df["Summary"] = results_df["PMCID"].map(summaries)
    progress_bar["value"] = 100
    window.update()

    # Generate the final PDF report and conclude the HERMES workflow
    print("[STATUS] Generating final PDF report...")
    progress_label2["text"] = "Generating final PDF report..."
    window.update()

    with open(log_file, 'a') as f:
        f.write("[STATUS] Generating final PDF report...\n")
    
    try:
        generate_pdf_report(date, terms, parsed, hits, keywords, hermes_mode, request, figures_dir, results_df, report_dir, log_file, idf)
        submit_label["text"] = ""
        progress_label2["text"] = "PDF report saved successfully."
        print("[INFO] PDF report saved successfully.")
    except:
        progress_label2["text"] = "[ERROR] Could not generate PDF report."
        with open(log_file, 'a') as f:
            f.write("[ERROR] Could not generate PDF report.\n")
    
    progress_bar["value"] = 0

    # Report total wall-clock runtime for the complete HERMES execution
    elapsed_total = time.perf_counter() - t0
    print("[INFO] Total run time: {elapsed_total:.2f} s".format(elapsed_total=elapsed_total))

    # Update peak memory usage
    update_peak_memory()
    print(f"[INFO] Peak memory usage: {peak_rss / 1024**2:.1f} MiB")

    # Update log
    with open(log_file, 'a') as f:
        f.write("[INFO] Total run time: {elapsed_total:.2f} s\n".format(elapsed_total=elapsed_total))
        f.write(f"[INFO] Peak memory usage: {peak_rss / 1024**2:.1f} MiB\n")

    # Exit pipeline
    print("[STATUS] Exiting HERMES.")
    window.update() 


##########################################
# Tkinter Graphical User Interface (GUI) #
##########################################

# Initialize and configure the main Tkinter window and application-wide GUI styling
window = tk.Tk()
window.minsize(1000, 650)
window.configure(bg='white')
window.title("Hermes (v1.2) - Open Source Literature Mining")
style = ttk.Style(window)
style.theme_use('default')
style.configure("custom.Horizontal.TProgressbar", troughcolor='white', background='#1F4063', thickness=20)
default_font = tkFont.nametofont("TkDefaultFont")
default_font.configure(family="Arial", size=9)

# Create the main container frame for GUI layout
frame = tk.Frame(window, bg="white")
frame.pack(padx=10, pady=10, anchor="nw")

# Header section: application branding and introductory description
img_logo = tk.PhotoImage(file=logo_path)
logo = tk.Label(master = frame, image=img_logo, bg='white', fg='white')
logo.grid(row=0, column=0, sticky="w")
greeting = tk.Label(master = frame, text = "Hermes is an open-source mining tool for open-access literature, enabling scoring, and ranking of full-text articles based on customizable keyword sets and other relevance metrics.", bg='white', fg='black')
greeting.grid(row=2, column=0, sticky="w")

# Request title input: define and validate a unique identifier for report generation
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

# Entrez server identification: email input required for PubMed Central queries
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

# PubMed Central query definition: search terms and optional PubMed filters
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

# Keyword configuration: define scoring keywords and optional strict mode
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

# Report size configuration: define the number of articles to include in the final report
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

# Submission controls
submit_frame = tk.Frame(frame, bg="white")
submit_frame.grid(row=19, column=0, sticky="w", pady=(15,0))
submit_button = tk.Button(master = submit_frame, text = "Submit Query", command = submit_query, bg='white', fg="black", relief='raised')
submit_button.grid(row=1, column=0, sticky="w", padx=2)
submit_label = tk.Label(master = submit_frame, text = "", bg='white', fg="black")
submit_label.grid(row=1, column=1, sticky="w", padx=3)

# Execution status and progress indicators
progress_bar = ttk.Progressbar(master = frame, orient = tk.HORIZONTAL, style="custom.Horizontal.TProgressbar", mode = "determinate", length = 1000)
progress_bar.grid(row=20, column=0, sticky="w", padx=2, pady=10)
progress_label0 = tk.Label(master = frame, text = "", bg='white', fg='black')
progress_label0.grid(row=21, column=0, sticky="w")
progress_label = tk.Label(master = frame, text = "", bg='white', fg='black')
progress_label.grid(row=22, column=0, sticky="w")
progress_label2 = tk.Label(master = frame, text = "", bg='white', fg='black')
progress_label2.grid(row=23, column=0, sticky="w")

# Output directory selection
my_dir = tk.StringVar(master = window, value = "")

# Enter the Tkinter event loop and wait for user interaction
window.mainloop()