#!/usr/bin/env python3

####################################################################
# Hermes v1.2 - Open-source mining tool for open-access literature #
# 2025-12-18                                                       #
# Written by Julien Charest & Katarina Priselac                    #
####################################################################

from os import path, mkdir
import sys


def ressource_path(relative_path):

    """Resolve the path to application resource."""

    try:
        base_path = sys._MEIPASS
    except AttributeError:
        base_path = path.abspath(".")
    return path.join(base_path, relative_path)


def create_reports_folder(working_directory, request_id, date):

    """Create a structured output directory for a HERMES run, including subfolders for figures and results."""

    report_dir = path.join("{working_directory}".format(working_directory = working_directory), "{date}_{request_id}".format(request_id = request_id, date = date))
    figures_dir = path.join("{report_dir}".format(report_dir = report_dir), "figures")
    results_dir = path.join("{report_dir}".format(report_dir = report_dir), "results")

    if not path.exists(report_dir):
        mkdir(report_dir)
    
    if not path.exists(figures_dir):
        mkdir(figures_dir)

    if not path.exists(results_dir):
        mkdir(results_dir)


def filter_string(string):

    """Remove unsupported characters from a string to ensure safe output in logs, files, and reports."""

    new_string = ""
    for i in string:
        if i in "abcdefghijklmnopqrstuvwxyz" + "abcdefghijklmnopqrstuvwxyz".upper() + " 0123456789.:,;.()!?=-+/%":
            new_string += i
    return str(new_string)