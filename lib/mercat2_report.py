# -*- coding: utf-8 -*-
"""mercat2_report.py: Module for creating tables, figures, and html reports
"""

import base64
import pkg_resources as pkg
import time
import numpy as np
import pandas as pd
import plotly.graph_objs as go
import plotly.express as px
from sklearn.decomposition import PCA
from functools import reduce
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
def PCA_plot(dfPCA):
    '''Creates a 3D plotly scatter plot of the PCA of the given Dataframe.

    Parameters:
        dfPCA (Dataframe): A pandas dataframe with samples and kmer counts

    Returns:
        plotly fig: A plotly 3D scatter plot of the PCA of the dataframe.
    '''

    df_merged2 = dfPCA
    result = df_merged2.replace(np.nan, 0)
    pivoted = result.T
    res = pivoted.rename(columns=pivoted.iloc[0])
    res1 = res.drop(res.index[0])
    pca = PCA(n_components=3, svd_solver='randomized')
    X_train = pca.fit_transform(res1)
    labels = {
    str(i): f"PC {i+1} ({var:.1f}%)"
    for i, var in enumerate(pca.explained_variance_ratio_ * 100)
    }
    figPCA = px.scatter_3d(
        X_train, x=0, y=1, z=2, color=res1.index,
        labels=labels,
        template="plotly_white"
    )
    figPCA.update_layout(font=dict(color="Black"),
        margin=dict(l=0, r=0, t=0, b=0),
        )
    
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
