# -*- coding: utf-8 -*-
"""mercat2_figures.py: Module for creating tables, figures
"""

import os
from pathlib import Path
import psutil
import re
import base64
import gzip
import pkg_resources as pkg
import pandas as pd
import plotly.graph_objs as go
import plotly.express as px
from sklearn.decomposition import PCA
from sklearn.decomposition import IncrementalPCA as iPCA
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


# RAM Usage
def mem_use():
    return round(psutil.virtual_memory().used/1024.0**3, 2)


# k-mer summary
def kmer_summary(tsv_file: dict):
    '''Creates a plotly bar graph with the k-mer counts of the top k-mers among all samples.

    Parameters: ***TODO: UPDATE THIS***
        df_all_samples (dict[name: Dataframe]): A dictionary with the Dataframes of samples containing the counts of k-mers.
        
    Returns:
        tuple: A tuple with two plotly figures: a barplot of the top k-mers and a plotly table figure labeling the k-mers.
    '''

    numKmers = 5

    top_kmers = list()
    def line_avg(str_list):
        return sum([int(x) for x in str_list[1:]]) / (len(str_list)-1)

    with open(tsv_file, 'r') as reader:
        header = reader.readline().strip().split(sep='\t')
        for line in reader:
            line = line.strip().split(sep='\t')
            if len(top_kmers) < numKmers:
                top_kmers.append(line)
            else:
                top_kmers.sort(key = lambda x: line_avg(x), reverse=False) # sort in ascending order by first element
                if line_avg(line) > line_avg(top_kmers[0]):
                    top_kmers[0] = line

    df = pd.DataFrame(top_kmers, columns=header).set_index('k-mer', drop=True).astype(int)
    df = df.T.reset_index().rename(columns=dict(index='sample'))
    df = df.melt(id_vars=['sample'], var_name='k-mer', value_name='count')#.sort_values(by='count', ascending=False)
    df['label'] = pd.Categorical(df['k-mer']).codes + 1
    df['label'] = df['label'].apply(lambda x: f'k-mer-{x}')
    df.sort_values(by=['label','count'], inplace=True, ascending=[True,False])

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

    return (fig, figTable)


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
def plot_sample_metrics(protein_samples: dict, tsv_out):
    '''Creates a plotly bar graph of the protein metrics.

    Parameters:
        protein_list (dict[str: list[str]]): A dictionary with lists of sequence files

    Returns:
        plotly fig: A plotly barplot figure of the protein metrics from the samples dictionary.
    '''

    with open(tsv_out, 'w') as writer:
        print('Sample', 'seq_name', 'length', 'PI', 'MW', 'Hydro', sep='\t', file=writer)

    figures = dict()
    for basename,files in protein_samples.items():
        for file in files:
            dfMetrics = pd.DataFrame()
            reader = gzip.open(file, 'rt') if Path(file).suffix=='.gz' else open(file, 'r')
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
                    if len(sequence):
                        dfMetrics.at[name, 'Name'] = name.split()[0]
                        dfMetrics.at[name, 'Length'] = len(sequence)
                        dfMetrics.at[name, 'PI'] = mercat2_metrics.predict_isoelectric_point_ProMoST(sequence)
                        dfMetrics.at[name, 'MW'] = mercat2_metrics.calculate_MW(sequence)
                        dfMetrics.at[name, 'Hydro'] = mercat2_metrics.calculate_hydro(sequence)
                    else:
                        print("WARNING: Empty Sequence:", basename, name, sequence)
                    continue #already got next line, next item in loop
                line = reader.readline()
            reader.close()

            dfMetrics.sort_values(by='Length', ascending=False, inplace=True)
            dfMetrics.to_csv(tsv_out, sep='\t', mode='a', index=True, header=False)

            figures[f"{basename}_PI"] = px.bar(dfMetrics, x='Length', y='PI', template="plotly_white",)
            figures[f"{basename}_PI"].update_layout(font=dict(color="Black"))

            figures[f"{basename}_MW"] = px.bar(dfMetrics, x='Length', y='MW', template="plotly_white",)
            figures[f"{basename}_MW"].update_layout(font=dict(color="Black"))

            figures[f"{basename}_Hydro"] = px.bar(dfMetrics, x='Length', y='Hydro', template="plotly_white",)
            figures[f"{basename}_Hydro"].update_layout(font=dict(color="Black"))

#            dfMetrics = dfMetrics.melt(id_vars=['Name','Length'], var_name='Metric', value_name='Value')
#            figures[f"{basename}_Metrics"] = px.bar(dfMetrics, x='Length', y='Value', facet_row='Metric', template="plotly_white")
#            figures[f"{basename}_Metrics"].update_yaxes(matches=None)
#            figures[f"{basename}_Metrics"].update_layout(font=dict(color="Black"))

    return figures


# PCA
def plot_PCA(tsv_file:str, out_path:str, lowmem=None, class_file=None, DEBUG=False):
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
    if DEBUG:
        print(f"\nTime to get name list: {round(timeit.default_timer() - start_time,2)} seconds")
        print(f"Virtual Memory {mem_use()}GB")

    if lowmem is None and len(names) > chunk_size:
        lowmem = True
    else:
        lowmem = False

    print("Using Incremental PCA:", lowmem)
    
    name_iter = iter(names)
    pca = iPCA(n_components=3, batch_size=100) if lowmem else PCA(n_components=3)
    XDF = pd.read_csv(tsv_file, sep='\t', index_col=0, chunksize=chunk_size) if lowmem else pd.read_csv(tsv_file, sep='\t', index_col=0)
    if DEBUG:
        print(f"\nTime to read TSV file: {round(timeit.default_timer() - start_time,2)} seconds")
        print(f"Virtual Memory {mem_use()}GB")
    
    if lowmem:
        # Incremental PCA
        # Partial Fit
        for chunk in XDF:
            pca.partial_fit(chunk)
        if DEBUG:
            print(f"\nTime for iPCA partial_fit: {round(timeit.default_timer() - start_time,2)} seconds")
            print(f"Virtual Memory {mem_use()}GB")

        # Transform
        XDF = pd.read_csv(tsv_file, sep='\t', index_col=0, chunksize=chunk_size)
        with open(pca_tsv, 'w') as pca_out:
            print('sample', 'PC1', 'PC2', 'PC3', sep='\t', file = pca_out)
            for chunk in XDF:
                for row in pca.transform(chunk):
                    pca_out.write(next(name_iter))
                    for c in row:
                        pca_out.write(f'\t{c}')
                    pca_out.write('\n')
        if DEBUG:
            print(f"\nTime for iPCA transform: {round(timeit.default_timer() - start_time,2)} seconds")
            print(f"Virtual Memory {mem_use()}GB")
    else:
        # Standard PCA
        XDF = pca.fit_transform(XDF)
        if DEBUG:
            print(f"\nTime to compute PCA fit_transform: {round(timeit.default_timer() - start_time,2)} seconds")
            print(f"Virtual Memory {mem_use()}GB")

        with open(pca_tsv, 'w') as pca_out:
            print('sample', 'PC1', 'PC2', 'PC3', sep='\t', file = pca_out)
            for row in XDF:
                pca_out.write(next(name_iter))
                for c in row:
                    pca_out.write(f'\t{c}')
                pca_out.write('\n')
        if DEBUG:
            print(f"Time to save PCA TSV file: {round(timeit.default_timer() - start_time,2)} seconds")
            print(f"Virtual Memory {mem_use()}GB")


    start_time = timeit.default_timer()
    XDF = pd.read_csv(pca_tsv, sep='\t', index_col=0)

    if class_file:
        df_tax = pd.read_csv(class_file, sep='\t', index_col=0, names=['class'])
        XDF['class'] = XDF.index.map(df_tax['class']).fillna('NA')
        color_col = 'class'
        if DEBUG:
            print(f"\nTime to load CLASS file: {round(timeit.default_timer() - start_time,2)} seconds")
            print(f"Virtual Memory {mem_use()}GB")
    else:
        color_col = names

    # Plotly PCA
    start_time = timeit.default_timer()
    labels_axis = {
        f'PC{i}' : f"PC {i} ({val:.1f}%)"
        for i,val in enumerate(pca.explained_variance_ratio_ * 100, start=1)
        }

    figPCA = px.scatter_3d(
        XDF, x='PC1', y='PC2', z='PC3', color=color_col,
        labels=labels_axis, template="plotly_white"
    )
    figPCA.update_layout(font=dict(color="Black"),
        margin=dict(l=0, r=0, t=0, b=0),)
    figPCA.update_layout(scene = dict(
        xaxis = dict(
            backgroundcolor="White",
            gridcolor="LightGray",
            showbackground=False,
            zerolinecolor="LightGray"
            ),
        yaxis = dict(
            backgroundcolor="White",
            gridcolor="LightGray",
            showbackground=False,
            zerolinecolor="LightGray"
            ),
        zaxis = dict(
            backgroundcolor="White",
            gridcolor="LightGray",
            showbackground=False,
            zerolinecolor="LightGray"
            ),
        ))

    figPCA.write_image(f"{out_path}/pca{'_incremental' if lowmem else ''}.png")
    
    print(f"Time to compute 3D PCA: {round(timeit.default_timer() - start_time,2)} seconds")
    if DEBUG:
        print(f"Virtual Memory {mem_use()}GB")

    start_time = timeit.default_timer()
    figPCA2d = None
    if pca.explained_variance_ratio_[2] * 100 < 1:
        figPCA2d = px.scatter(
            XDF, x='PC1', y='PC2', color=color_col,
            labels=labels_axis, template="plotly_white"
        )
        figPCA2d.update_layout(font=dict(color="Black"),
            margin=dict(l=0, r=0, t=0, b=0),)
        figPCA2d.write_image(f"{out_path}/pca2D{'_incremental' if lowmem else ''}.png")
        print(f"Time to compute 2D PCA: {round(timeit.default_timer() - start_time,2)} seconds")
        if DEBUG:
            print(f"Virtual Memory {mem_use()}GB")

    return figPCA, figPCA2d
