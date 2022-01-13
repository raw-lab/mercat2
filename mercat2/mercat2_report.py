import os
import base64
import pkg_resources as pkg
import time
import numpy as np
import pandas as pd
import itertools
import plotly.graph_objs as go
import matplotlib.pyplot as plt
import plotly.express as px
from sklearn.decomposition import PCA
from functools import reduce
import dominate
from dominate.tags import *
from dominate.util import raw

from mercat2 import mercat2_metrics


# Global Stylesheet
STYLESHEET = pkg.resource_stream('mercat2', 'data/style.css').read().decode()
LOGO = pkg.resource_stream('mercat2', 'data/logo.jpg').read()
LOGO = base64.b64encode(LOGO).decode()

# TODO: Option for how to include plotly.js.
# False uses script in <head>, 'cdn' loads from internet. # Can I use both???
PLOTLY_SOURCE = 'cdn'


# k-mer summary
def kmer_summary(df_all_samples):
    df_list = []
    for name,value in df_all_samples.items():
        df_list.append(value['Count'].rename(name))
    df_merged = pd.concat(df_list, axis=1, keys=[s.name for s in df_list]).fillna(0)
    df_merged['Total'] = df_merged.sum(axis=1)
    df_merged['Mean'] = df_merged[[s.name for s in df_list]].mean(axis='columns')
    df_merged.sort_values(by='Mean', axis='index', ascending=False, inplace=True)
    df = df_merged.iloc[0:2][[s.name for s in df_list]].T.reset_index().rename(columns=dict(index='sample'))
    df = df.melt(id_vars=['sample'], var_name='k-mer', value_name='count').sort_values(by='count', ascending=False)
    df['label'] = pd.Categorical(df['k-mer']).codes + 1
    df['label'] = df['label'].apply(lambda x: f'kmer-{x}')
    
    fig = px.bar(df, x='sample', y='count', text='count', color='sample', facet_row='label')
    fig.for_each_annotation(lambda a: a.update(text=a.text.replace('label=', '')))
    
    df = df[['k-mer', 'label']].drop_duplicates()

    # Add table to label kmers
    figTable = go.Figure(
        data=[go.Table(cells=dict(values=[df['label'].tolist(), df['k-mer'].tolist()]))],
        layout=go.Layout(margin=go.layout.Margin(
            l=0, #left margin
            r=0, #right margin
            b=0, #bottom margin
            t=0, #top margin
        ))
    )
    figTable.update_layout(height=100)
    return [fig, figTable]


# GC Plot kmers
def GC_plot_kmer(df_all_samples):
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
    fig = px.bar(df, x='k-mer', y='GC', text='GC')
    return fig


# GC Plot samples
def GC_plot_sample(gc_content: dict):
    print(gc_content)
    df = pd.DataFrame.from_dict(data=gc_content, orient='index', columns=['GC Content'])
    fig = px.bar(df)
    return fig


# Plot Protein Metrics samples
def Metric_plot_samples(kmer_list):
    df = pd.DataFrame(0.0, index=kmer_list, columns=['Count',"PI","MW","Hydro"])
    for k in kmer_list:
        df.at[k,'Count'] = k
        df.at[k,'PI'] = mercat2_metrics.predict_isoelectric_point_ProMoST(k)
        df.at[k,'MW'] = mercat2_metrics.calculate_MW(k)
        df.at[k,'Hydro'] = mercat2_metrics.calculate_hydro(k)
    print(df)
    fig = px.bar(df)
    return fig


# PCA
def PCA_plot(dfPCA):
    df_merged2 = dfPCA
    result = df_merged2.replace(np.nan, 0)
    pivoted = result.T
    res = pivoted.rename(columns=pivoted.iloc[0])
    res1 = res.drop(res.index[0])
    pca = PCA(n_components=3,svd_solver='randomized')
    X_train = pca.fit_transform(res1)
    labels = {
    str(i): f"PC {i+1} ({var:.1f}%)"
    for i, var in enumerate(pca.explained_variance_ratio_ * 100)
    }
    figPCA = px.scatter_3d(
        X_train, x=0, y=1, z=2, color=res1.index,
        labels=labels
    )
    figPCA.update_layout({
    'plot_bgcolor' : '#7f7f7f',
    'paper_bgcolor': '#FFFFFF',
    'paper_bgcolor': "rgba(0,0,0,0)",
    'plot_bgcolor' : '#7f7f7f',
    
    })
    figPCA.update_layout(scene = dict(
                        xaxis = dict(
                            backgroundcolor="rgb(255,255, 255)",
                            gridcolor="black",
                            showbackground=True,
                            zerolinecolor="black",),
                        yaxis = dict(
                            backgroundcolor="rgb(255,255, 255)",
                            gridcolor="black",
                            showbackground=True,
                            zerolinecolor="black"),
                        zaxis = dict(
                            backgroundcolor="rgb(255,255, 255)",
                            gridcolor="black",
                            showbackground=True,
                            zerolinecolor="black",),),
                    
                    )
    
    pca1 = PCA(n_components=5,svd_solver='randomized')
    X_train = pca1.fit_transform(res1)
    per_var = np.round(pca1.explained_variance_ratio_*100, decimals=1)
    labels = ['PC'+str(i) for i in range(1,len(per_var)+1)]
    fig = plt.bar(range(1,len(per_var)+1), per_var, tick_label=labels)
    plt.xlabel('Principal Components',fontsize=15)
    plt.ylabel('Variance %',fontsize=15)
    #plt.savefig(subdir+'/'+'plots'+'/'+'PCA_plot_ScreePlot'+".png")
    return figPCA


# Write HTML Report
def write_html(outfile, figPlots, tsv_stats):
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
