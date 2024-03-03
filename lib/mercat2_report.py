# -*- coding: utf-8 -*-
"""mercat2_report.py: Module for creating html reports
"""

import os
import resource
import psutil
import re
import base64
import pkg_resources as pkg
import time
import dominate
from dominate.tags import *
from dominate.util import raw

from mercat2_lib import mercat2_metrics


# Global Stylesheet
STYLESHEET = pkg.resource_stream('mercat2_lib', 'data/style.css').read().decode()
LOGO = pkg.resource_stream('mercat2_lib', 'data/logo.jpg').read()
LOGO = base64.b64encode(LOGO).decode()

# TODO: Option for how to include plotly.js.
# False uses script in <head>, 'cdn' loads from internet. # Can I use both???
PLOTLY_SOURCE = 'cdn'


# RAM Usage
def mem_use():
    return round(psutil.virtual_memory().used/1024.0**3, 2)


# Write HTML Report
def write_html(outfile:str, figPlots:dict, tsv_stats:dict):
    '''Creates an html report with the figures and tsv files.

    Parameters:
        outfile (str): A path to the html file to save the report.
        figPlots (dict): A dictionary of plotly figures.
        tsv_stats (dict): A dictionary of paths to TSV files to add to the downloads section.

    Returns:
        None
    '''

    # WRITE TO FORMATED HTML
    with dominate.document(title='K-Mer Report') as doc:
        # Header
        with doc.head:
            meta(charset="utf-8")
            script(type="text/javascript", src="plotly-2.0.0.min.js")
            with style(type="text/css"):
                raw('\n'+STYLESHEET)
        # Body of HTML
        with div(cls="document", id="mercat-summary"):
            # Title bar at the top
            with h1(cls="title"):
                img(src=f"data:image/png;base64,{LOGO}", height="40")
                a("MERCAT2", cls="reference external", href="https://github.com/raw-lab/mercat2")
                raw(" - Summary")
            # Table of contents - left panel
            with div(cls="contents topic", id="contents"):
                with ul(cls="simple"):
                    li(a("Summary", cls="reference internal", href="#summary"))
                    with ul():
                        for key in figPlots.keys():
                            li(a(f"{key}", cls="reference internal", href=f"#{key}"))
                    li(a("Downloads", cls="reference internal", href="#downloads"))
            # Main panel
            with div(h1("Summary"), cls="section", id="summary"):
                for key,figures in figPlots.items():
                    with div(h2(f"{key}"), cls="section", id=f"{key}"):
                        try:
                            raw(figures.to_html(full_html=False, include_plotlyjs=PLOTLY_SOURCE))
                        except:
                            for fig in figures:
                                raw(fig.to_html(full_html=False, include_plotlyjs=PLOTLY_SOURCE))
            # Downloads at the bottom
            with div(cls="section", id="downloads"):
                h1("Downloads")
                with div(cls="docutils container", id="attachments"):
                    with blockquote():
                        with div(cls="docutils container", id="table-1"):
                            with dl(cls="docutils"):
                                for key,value in tsv_stats.items():
                                    data_URI = f"data:text/tab-separated-values;base64,{value}"
                                    dt(key)
                                    dd(a(key+".tsv", href=data_URI, download=key+".tsv", draggable="true"))
                div(time.strftime("%Y-%m-%d", time.localtime()), cls="docutils container", id="metadata")

    with open(outfile, 'w') as writer:
        writer.write(doc.render())
    return


# Merge TSV Files
def merge_tsv(tsv_list:dict, out_file:os.PathLike):
    FLIMIT = len(tsv_list) >= resource.getrlimit(resource.RLIMIT_NOFILE)[0]
    names = sorted(list(tsv_list.keys()))
    file_list = dict()
    header = ""
    for name in names:
        if FLIMIT:
            with open(tsv_list[name]) as f:
                head = f.readline() # skip header
                if not header:
                    header = head.split('\t')[0]
                file_list[name] = f.tell()
        else:
            file_list[name] = open(tsv_list[name])
            head = file_list[name].readline() # skip header
            if not header:
                header = head.split('\t')[0]
    with open(out_file, 'w') as writer:
        print(header, '\t'.join(names), sep='\t', file=writer)
        lines = dict()
        kmers = set()
        for name in names:
            if FLIMIT:
                with open(tsv_list[name]) as f:
                    f.seek(file_list[name])
                    lines[name] = f.readline().split()
                    file_list[name] = f.tell()
            else:
                lines[name] = file_list[name].readline().split()
            kmers.add(lines[name][0])
        kmer = sorted(kmers)[0]
        while True:
            line = [kmer]
            kmers = set()
            for name in names:
                if not lines[name]:
                    line.append('0')
                elif lines[name][0] > kmer:
                    line.append('0')
                else:
                    line.append(lines[name][1])
                    if FLIMIT:
                        with open(tsv_list[name]) as f:
                            f.seek(file_list[name])
                            lines[name] = f.readline().split()
                            file_list[name] = f.tell()
                    else:
                        lines[name] = file_list[name].readline().strip('\n').split('\t')
                        lines[name] = [x for x in lines[name] if len(x) > 0]
                    if lines[name]:
                        kmers.add(lines[name][0])
            print('\t'.join(line), file=writer)
            if not kmers:
                break
            kmer = sorted(kmers)[0]
    if FLIMIT:
        for name in names:
            file_list[name].close()
    return


# Merge TSV Files, Transposed
def merge_tsv_T(tsv_list:dict, out_file:os.PathLike):
    names = sorted(list(tsv_list.keys()))#, key=lambda x: x.lower())
    header = set()
    #file_list = dict()
    for name in names:
        with open(tsv_list[name]) as reader:
            #file_list[name] = open(tsv_list[name], 'r')
            #file_list[name].readline() # skip header
            reader.readline() # skip header
            for line in reader:
                kmer = line.split()[0]
                header.add(kmer)
    
    header = list(header)
    with open(out_file, 'w') as writer:
        print('sample', '\t'.join(header), sep='\t', file=writer)
        for name in names:
            with open(tsv_list[name]) as reader:
                #file_list[name].seek(0)
                #file_list[name].readline() # skip header
                reader.readline() # skip header
                counts = dict()
                #for line in file_list[name]:
                for line in reader:
                    kmer,count = line.split()
                    counts[kmer] = count
                #file_list[name].close()
            writer.write(name)
            for kmer in header:
                if kmer in counts:
                    writer.write(f'\t{counts[kmer]}')
                else:
                    writer.write(f'\t0')
            writer.write('\n')
    return
