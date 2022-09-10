# -*- coding: utf-8 -*-
"""mercat2_report.py: Module for creating tables, figures, and html reports
"""

import os
import re
import base64
import pkg_resources as pkg
import time
import pandas as pd
import plotly.graph_objs as go
import plotly.express as px
from sklearn.decomposition import PCA
from sklearn.decomposition import IncrementalPCA as iPCA
import dominate
from dominate.tags import *
from dominate.util import raw
import timeit

from mercat2_lib import mercat2_metrics


# Global Stylesheet
STYLESHEET = pkg.resource_stream('mercat2_lib', 'data/style.css').read().decode()
LOGO = pkg.resource_stream('mercat2_lib', 'data/logo.jpg').read()
LOGO = base64.b64encode(LOGO).decode()

# TODO: Option for how to include plotly.js.
# False uses script in <head>, 'cdn' loads from internet. # Can I use both???
PLOTLY_SOURCE = 'cdn'


# k-mer summary
def kmer_summary(df_all_samples: dict):
    '''Creates a plotly bar graph with the k-mer counts of the top 2 k-mers among all samples.

    Parameters:
        df_all_samples (dict[name: Dataframe]): A dictionary with the Dataframes of samples containing the counts of k-mers.

    Returns:
        tuple: A tuple with two plotly figures: a barplot of the top k-mers and a plotly table figure labeling the k-mers.
    '''

    numKmers = 5

    df_list = []
    for name,value in df_all_samples.items():
        df_list.append(value['Count'].rename(name))
    df_merged = pd.concat(df_list, axis=1, keys=[s.name for s in df_list]).fillna(0)
    df_merged['Total'] = df_merged.sum(axis=1)
    df_merged['Mean'] = df_merged[[s.name for s in df_list]].mean(axis='columns')
    df_merged.sort_values(by='Mean', axis='index', ascending=False, inplace=True)
    df = df_merged.iloc[0:numKmers][[s.name for s in df_list]].T.reset_index().rename(columns=dict(index='sample'))
    df = df.melt(id_vars=['sample'], var_name='k-mer', value_name='count').sort_values(by='count', ascending=False)
    df['label'] = pd.Categorical(df['k-mer']).codes + 1
    df['label'] = df['label'].apply(lambda x: f'kmer-{x}')
    
    fig = px.bar(df, x='sample', y='count', text='count', color='sample', facet_row='label', template="plotly_white")
    fig.for_each_annotation(lambda a: a.update(text=a.text.replace('label=', '')))
    fig.update_layout(font=dict(color="Black"))
    
    df = df[['k-mer', 'label']].drop_duplicates()

    # Add table to label kmers
    figTable = go.Figure(
        data=[go.Table(cells=dict(values=[df['label'].tolist(), df['k-mer'].tolist()]))],
        layout=go.Layout(margin=go.layout.Margin(l=0, r=0, b=0, t=0,
        ))
    )
    figTable.update_layout(height=100, font=dict(color="Black"), template='plotly_white')
    return [fig, figTable]


# GC Plot kmers
def GC_plot_kmer(df_all_samples:dict):
    '''Creates a plotly bar graph with the GC content of k-mers.

    Parameters:
        df_all_samples (dict[name: Dataframe]): A dictionary with the Dataframes of samples containing the counts of k-mers.

    Returns:
        plotly fig: A plotly barplot figure of the GC% Content of the top k-mers.
    '''

    #print("\nGC PLOT\n")
    df_list = []
    for name,value in df_all_samples.items():
        df_list.append(value['Count'].rename(name))
    df_merged = pd.concat(df_list, axis=1, keys=[s.name for s in df_list]).fillna(0)
    df_merged['Total'] = df_merged.sum(axis=1)
    df_merged['Mean'] = df_merged[[s.name for s in df_list]].mean(axis='columns')
    df_merged.sort_values(by='Mean', axis='index', ascending=False, inplace=True)
    #print(df_merged)
    df = df_merged.reset_index().rename(columns=dict(index='k-mer'))['k-mer'].iloc[0:2]
    def calc_gc(kmer:str):
        return 100.0 * (kmer.count('G') + kmer.count('C')) / len(kmer)
    df = pd.concat([df, df.apply(calc_gc)], axis=1, keys=['k-mer', 'GC'])
    #print("GC CONTENT\n", df)
    fig = px.bar(df, x='k-mer', y='GC', text='GC', template="plotly_white")
    fig.update_layout(font=dict(color="Black"))
    return fig


# GC Plot samples
def GC_plot_sample(gc_content: dict):
    '''Creates a plotly bar graph of the GC content.

    Parameters:
        gc_content (dict[str: float]): A dictionary with GC content data.

    Returns:
        plotly fig: A plotly barplot figure of the GC% Content from the gc_content dictionary.
    '''

    df = pd.DataFrame.from_dict(data=gc_content, orient='index', columns=['GC Content'])
    fig = px.bar(df, template="plotly_white",
        labels={'index':'Sample', 'value':'GC percent'})
    fig.update_layout(font=dict(color="Black"))
    return fig


# Protein Metrics Plot samples
def plot_sample_metrics(protein_samples: dict):
    '''Creates a plotly bar graph of the protein metrics.

    Parameters:
        protein_list (dict[str: list[str]]): A dictionary with lists of sequence files

    Returns:
        plotly fig: A plotly barplot figure of the protein metrics from the samples dictionary.
    '''

    dfMetrics = pd.DataFrame()
    for basename,files in protein_samples.items():
        prot_count = 0
        prot_PI = 0.0
        prot_MW = 0.0
        prot_Hydro = 0.0
        for file in files:
            with open(file) as reader:
                line = reader.readline()
                while line:
                    line = line.strip()
                    line = line.rstrip('*')
                    if line.startswith('>'):
                        name = line[1:]
                        sequence = ""
                        line = reader.readline()
                        while line:
                            line = line.strip()
                            line = line.rstrip('*')
                            if line.startswith('>'):
                                break
                            sequence += line
                            line = reader.readline()
                        prot_count += 1
                        prot_PI = mercat2_metrics.predict_isoelectric_point_ProMoST(sequence)
                        prot_MW = mercat2_metrics.calculate_MW(sequence)
                        prot_Hydro = mercat2_metrics.calculate_hydro(sequence)
                        continue #already got next line, next item in loop
                    line = reader.readline()
        dfMetrics.at[basename, 'PI Avg'] = prot_PI / prot_count
        dfMetrics.at[basename, 'MW Avg'] = prot_MW / prot_count
        dfMetrics.at[basename, 'Hydro Avg'] = prot_Hydro / prot_count

    figPI = px.bar(dfMetrics, x=dfMetrics.index, y='PI Avg', template="plotly_white",
        labels={'index':'Sample', 'PI Avg':'PI Average'})
    figPI.update_layout(font=dict(color="Black"))
    
    figMW = px.bar(dfMetrics, x=dfMetrics.index, y='MW Avg', template="plotly_white",
        labels={'index':'Sample', 'MW Avg':'MW Average'})
    figMW.update_layout(font=dict(color="Black"))
    
    figHydro = px.bar(dfMetrics, x=dfMetrics.index, y='Hydro Avg', template="plotly_white",
        labels={'index':'Sample', 'Hydro Avg':'Hydropathy Average'})
    figHydro.update_layout(font=dict(color="Black"))
    return ([figPI, figMW, figHydro], dfMetrics)


# PCA
def plot_PCA(tsv_file:str, out_path:str, lowmem=None):
    '''Creates a 3D plotly scatter plot of the PCA of the given TSV File.

    Parameters:
        tsv_file (str): path to a TSV file with samples and kmer counts
        out_path (str): path to folder to save the PCA image and transformed PCA components

    Returns:
        plotly fig: A plotly 3D scatter plot of the PCA of the dataframe.
    '''

    start_time = timeit.default_timer()

    pca_tsv = os.path.join(out_path, 'pca.tsv')
    chunk_size = 1000
    names = list()
    with open(tsv_file) as reader:
        reader.readline() # skip header
        for line in reader:
            names.append(re.sub(r'_protein', '', line.split()[0]))
    print(f"Time to get name list: {round(timeit.default_timer() - start_time,2)} seconds")

    if lowmem is None and len(names) > chunk_size:
        lowmem = True
    else:
        lowmem = False

    print("Using low memory:", lowmem)
    
    name_iter = iter(names)
    pca = iPCA(n_components=3, batch_size=100) if lowmem else PCA(n_components=3)
    df = pd.read_csv(tsv_file, sep='\t', index_col=0, chunksize=chunk_size) if lowmem else pd.read_csv(tsv_file, sep='\t', index_col=0)
    print(f"Time to read TSV file: {round(timeit.default_timer() - start_time,2)} seconds")
    
    if lowmem:
        # Incremental PCA
        # Partial Fit
        for chunk in df:
            pca.partial_fit(chunk)
        print(f"Time for iPCA partial_fit: {round(timeit.default_timer() - start_time,2)} seconds")
        # Transform
        df = pd.read_csv(tsv_file, sep='\t', index_col=0, chunksize=chunk_size)
        with open(pca_tsv, 'w') as pca_out:
            print('sample', 'PC1', 'PC2', 'PC3', sep='\t', file = pca_out)
            for chunk in df:
                for row in pca.transform(chunk):
                    pca_out.write(next(name_iter))
                    for c in row:
                        pca_out.write(f'\t{c}')
                    pca_out.write('\n')
        print(f"Time for iPCA transform: {round(timeit.default_timer() - start_time,2)} seconds")
    else:
        # Standard PCA
        X = pca.fit_transform(df)
        print(f"Time to compute PCA fit_transform: {round(timeit.default_timer() - start_time,2)} seconds")
        with open(pca_tsv, 'w') as pca_out:
            print('sample', 'PC1', 'PC2', 'PC3', sep='\t', file = pca_out)
            for row in X:
                pca_out.write(next(name_iter))
                for c in row:
                    pca_out.write(f'\t{c}')
                pca_out.write('\n')
        print(f"Time to save PCA TSV file: {round(timeit.default_timer() - start_time,2)} seconds")

    Xdf = pd.read_csv(pca_tsv, sep='\t', index_col=0)
    print(f"Time to read PCA TSV file: {round(timeit.default_timer() - start_time,2)} seconds")

    tax_file = '/projects/raw_lab/databases/GTDB/all_taxonomy_r207.tsv'
    tax_rank = ['Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
    df_tax = pd.read_csv(tax_file, sep='\t', index_col=0, names=['taxonomy'])
    #RS_GCF_000566285.1     d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli
    Xdf['taxonomy'] = Xdf.index.map(df_tax['taxonomy']).fillna('NA;NA;NA;NA;NA;NA;NA')
    Xdf[tax_rank] = Xdf['taxonomy'].str.split(';', expand=True)
    print(Xdf)

    for tax in tax_rank:
        Xdf[tax] = Xdf[tax].apply(lambda x: x[3:])
        print(f"Unique {tax}:", Xdf[tax].nunique())

    # Plotly PCA
    print("Creating Figures")
    labels_axis = {
        f'PC{i}' : f"PC {i} ({val:.1f}%)"
        for i,val in enumerate(pca.explained_variance_ratio_ * 100, start=1)
        }

    figPCA = px.scatter_3d(
        Xdf, x='PC1', y='PC2', z='PC3', color='Order',
        labels=labels_axis,
        template="plotly_white"
    )
    figPCA.update_layout(font=dict(color="Black"),
        margin=dict(l=0, r=0, t=0, b=0),)
    figPCA.write_json(f"{out_path}/pca{'_incremental' if lowmem else ''}.json")
    figPCA.write_image(f"{out_path}/pca{'_incremental' if lowmem else ''}.png")
    print(f"Time to compute Figures: {round(timeit.default_timer() - start_time,2)} seconds")
    
    if pca.explained_variance_ratio_[2] * 100 < 5:
        figPCA = px.scatter(
            Xdf, x='PC1', y='PC2', color='Order',
            labels=labels_axis,
            template="plotly_white"
        )
        figPCA.update_layout(font=dict(color="Black"),
            margin=dict(l=0, r=0, t=0, b=0),)
        figPCA.write_json(f"{out_path}/pca2D{'_incremental' if lowmem else ''}.json")
        figPCA.write_image(f"{out_path}/pca2D{'_incremental' if lowmem else ''}.png")
        print(f"Time to compute 2D Figures: {round(timeit.default_timer() - start_time,2)} seconds")

    return figPCA


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
                        if type(figures) is list:
                            for fig in figures:
                                raw(fig.to_html(full_html=False, include_plotlyjs=PLOTLY_SOURCE))
                        else:
                            raw(figures.to_html(full_html=False, include_plotlyjs=PLOTLY_SOURCE))
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
    names = sorted(list(tsv_list.keys()))
    file_list = dict()
    for name in names:
        file_list[name] = open(tsv_list[name])
        file_list[name].readline() # skip header
    with open(out_file, 'w') as writer:
        print("kmer", '\t'.join(names), sep='\t', file=writer)
        lines = dict()
        kmers = set()
        for name in names:
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
                    lines[name] = file_list[name].readline().split()
                    if lines[name]:
                        kmers.add(lines[name][0])
            print('\t'.join(line), file=writer)
            if not kmers:
                break
            kmer = sorted(kmers)[0]
    for name in names:
        file_list[name].close()
    return


# Merge TSV Files, Transposed
def merge_tsv_T(tsv_list:dict, out_file:os.PathLike):
    names = sorted(list(tsv_list.keys()))
    header = set()
    file_list = dict()
    for name in names:
        file_list[name] = open(tsv_list[name], 'r')
        file_list[name].readline() # skip header
        for line in file_list[name]:
            kmer = line.split()[0]
            header.add(kmer)
    
    header = list(header)
    with open(out_file, 'w') as writer:
        print('sample', '\t'.join(header), sep='\t', file=writer)
        for name in names:
            file_list[name].seek(0)
            file_list[name].readline() # skip header
            counts = dict()
            for line in file_list[name]:
                kmer,count = line.split()
                counts[kmer] = count
            file_list[name].close()
            writer.write(name)
            for kmer in header:
                if kmer in counts:
                    writer.write(f'\t{counts[kmer]}')
                else:
                    writer.write(f'\t0')
            writer.write('\n')
    return
