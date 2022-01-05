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


# Global Stylesheet
STYLESHEET = pkg.resource_stream('mercat2', 'data/style.css').read().decode()
LOGO = pkg.resource_stream('mercat2', 'data/logo.jpg').read()
LOGO = base64.b64encode(LOGO).decode()

# TODO: Option for how to include plotly.js.
# False uses script in <head>, 'cdn' loads from internet. # Can I use both???
PLOTLY_SOURCE = 'cdn'


# Scatter Plots
def scatter_plots(bif,xlab,res_df,kmerstring):

    axis_title_font_size = 20
    axis_tick_label_size = 18
    legend_font_size = 14
    marker_size = 10

    save_img = 0

    trace1 = go.Scatter(
        x=res_df[xlab].values,
        y=res_df.index.values,
        mode='markers'
    )

    data = go.Figure(trace1)
    layout = go.Layout(
        legend=dict(
            font=dict(
                # family='sans-serif',
                size=legend_font_size,
                color='black'
            ),
            # borderwidth=2
        ),
        autosize=True,
        height=632,
        # title='Baseline',
        width=1274,
        xaxis=dict(
            # autorange=True,
            fixedrange=False,
            title=xlab,
            # type='linear',
            # showgrid=False,
            ticks='inside',
            ticklen=8,
            tickwidth=2,
            tickcolor='#000',
            #tickvals=list(reversed(cores)),
            titlefont=dict(
                # family='Courier New, monospace',
                size=axis_title_font_size,
                color='black'
            ),
            tickfont=dict(
                # family='Old Standard TT, serif',
                size=axis_tick_label_size,
                color='black'
            ),
            showgrid=True,
            showline=True,
            mirror=True,
            zeroline=False,
            # gridcolor='black',
            # linecolor='black',
            gridwidth=2
        ),
        yaxis=dict(
            #type='log',
            autorange=True,
            fixedrange=False,
            ticks='inside',
            ticklen=8,
            tickwidth=2,
            tickcolor='#000',
            titlefont=dict(
                # family='Courier New, monospace',
                size=axis_title_font_size,
                color='black'
                # color='#7f7f7f'
            ),
            rangemode='normal',
            tickmode='linear',
            # tickwidth=4,
            tickfont=dict(
                # family='Old Standard TT, serif',
                size=axis_tick_label_size,
                color='black'
            ),
            showgrid=True,
            showline=True,
            mirror=True,
            zeroline=False,
            # gridcolor='black',
            # linecolor='black',
            gridwidth=2,
            title=kmerstring
        )
    )

    fig = go.Figure(data=data, layout=layout)
    #plot(fig, filename=bif + "_"+xlab+".html", auto_open=False)
    return fig


# Stacked Bar Plots
def stackedbar_plots(top10_all_samples, kmerstring):
    axis_title_font_size = 20
    axis_tick_label_size = 18
    legend_font_size = 14
    marker_size = 10

    topk10 = 10
    all_samples = list(top10_all_samples.keys())
    btraces = []
    kmer_percent = dict()

    kmernames = dict()

    for bif in top10_all_samples:
        res_df, total_freq_count = top10_all_samples[bif]
        index_vals = res_df.index.values
        kmernames[bif]=index_vals
        #topk10 = min(len(index_vals), 10)
        #If a small sample has less than 10kmers, ignore it
        if len(index_vals) < 10: continue

        kmer_percent[bif] = []
        for i in range(0, topk10):
            fr = res_df.loc[index_vals[i], 'Count']
            fr = round(100.0 * float(fr) / total_freq_count, 1)
            kmer_percent[bif].append(fr)

    ylist = dict()
    all_kmer_percents = (list(kmer_percent.values()))
    ylist[0] = list(itertools.repeat(0,len(all_samples)))
    for i in range(1,topk10+1):
        ylist[i] = [x[i-1] for x in all_kmer_percents]

    all_annotations = []

    for s in ylist:
        if s==0: continue
        trace1 = go.Bar(
            x=all_samples,
            y=ylist[s],
            name="kmer" + str(s)
        )
        btraces.append(trace1)

        annotations = []
        for index, (xi, yi) in enumerate(zip(all_samples, ylist[s])):
            prevy = 0
            for z in range(0, s):
                prevy += ylist[z][index]

            annotations.append(
                 dict(x=xi, y=yi+prevy,
                 text=str(kmernames[xi][s-1]+": " + str(yi)+"%"),
                 xanchor='center',
                 yanchor='bottom',
                 showarrow=False,
                 )
            )

        all_annotations.append(annotations)


    data = go.Figure(btraces)
    layout = go.Layout(
        barmode='stack',
        annotations=list(itertools.chain.from_iterable(all_annotations)),
        legend=dict(
            font=dict(
                # family='sans-serif',
                size=legend_font_size,
                color='black'
            ),
            # borderwidth=2
        ),
        autosize=True,
        height=632,
        # title='Baseline',
        width=1274,
        xaxis=dict(
            # autorange=True,
            fixedrange=False,
            title="Samples",
            #type='linear',
            # showgrid=False,
            ticks='inside',
            ticklen=8,
            tickwidth=2,
            tickcolor='#000',
            #tickvals=list(reversed(cores)),
            titlefont=dict(
                # family='Courier New, monospace',
                size=axis_title_font_size,
                color='black'
            ),
            tickfont=dict(
                # family='Old Standard TT, serif',
                size=axis_tick_label_size,
                color='black'
            ),
            showgrid=True,
            showline=True,
            mirror=True,
            zeroline=False,
            # gridcolor='black',
            # linecolor='black',
            gridwidth=2
        ),
        yaxis=dict(
            #type='linear',
            autorange=True,
            fixedrange=False,
            ticks='inside',
            ticklen=8,
            tickwidth=2,
            tickcolor='#000',
            titlefont=dict(
                # family='Courier New, monospace',
                size=axis_title_font_size,
                color='black'
                # color='#7f7f7f'
            ),
            rangemode='normal',
            tickmode='linear',
            # tickwidth=4,
            tickfont=dict(
                # family='Old Standard TT, serif',
                size=axis_tick_label_size,
                color='black'
            ),
            showgrid=True,
            showline=True,
            mirror=True,
            zeroline=False,
            # gridcolor='black',
            # linecolor='black',
            gridwidth=2,
            title=kmerstring
        )
    )

    fig = go.Figure(data=data, layout=layout)
    #plot(fig, filename=inp_folder + "_barchart_"+xlab+".html", auto_open=False)
    return fig


# PCA
def PCA_plot(subdir):


    parent_dir = subdir
    subject_dirs = [os.path.join(parent_dir, dir) for dir in os.listdir(parent_dir) if os.path.isdir(os.path.join(parent_dir, dir))]

    filelist = []
    for dir in subject_dirs:
        csv_files = [os.path.join(dir, csv) for csv in os.listdir(dir) if os.path.isfile(os.path.join(dir, csv)) and csv.endswith('summary.csv')]
        for file in csv_files:
            a=pd.read_csv(file)
            base=os.path.basename(file)
            y=os.path.splitext(base)[0]
            y = y.replace('_nucleotide_summary','').replace('_protein_summary','')
            a = a[['k-mers','Count']]
            a = a.rename(columns={'Count':y })
            filelist.append(a)
            

    df_merged2 = reduce(lambda  left,right: pd.merge(left,right,on=['k-mers'],
                                                how='outer'), filelist)
    result=df_merged2.replace(np.nan, 0)
    pivoted=result.T
    res=pivoted.rename(columns=pivoted.iloc[0])
    res1=res.drop(res.index[0])
    pca = PCA(n_components=3,svd_solver='randomized')
    X_train= pca.fit_transform(res1)
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
    
    #plot(figPCA, filename=subdir+'/plots/'+'PCA_plot'+".html", auto_open=False)
    pca1 = PCA(n_components=5,svd_solver='randomized')
    X_train = pca1.fit_transform(res1)
    per_var = np.round(pca1.explained_variance_ratio_*100, decimals=1)
    labels = ['PC'+str(i) for i in range(1,len(per_var)+1)]
    fig = plt.bar(range(1,len(per_var)+1), per_var, tick_label=labels)
    plt.xlabel('Principal Components',fontsize=15)
    plt.ylabel('Variance %',fontsize=15)
    plt.savefig(subdir+'/'+'plots'+'/'+'PCA_plot_ScreePlot'+".png")
    return figPCA


# Write HTML Report
def write_html(outfile, figPlots, tsv_stats):
    # WRITE TO FORMATED HTML
    with dominate.document(title='K-Mer Report') as doc:
        with doc.head:
            meta(charset="utf-8")
            script(type="text/javascript", src="plotly-2.0.0.min.js")
            with style(type="text/css"):
                raw('\n'+STYLESHEET)
        with div(cls="document", id="mercat-summary"):
            with h1(cls="title"):
                img(src=f"data:image/png;base64,{LOGO}", height="40")
                a("MERCAT2", cls="reference external", href="https://github.com/raw-lab/mercat2")
                raw(" - Statistical Summary")
            with div(cls="contents topic", id="contents"):
                with ul(cls="simple"):
                    li(a("Summary", cls="reference internal", href="#summary"))
                    with ul():
                        for key in figPlots.keys():
                            li(a(f"{key}", cls="reference internal", href=f"#{key}"))
                    li(a("Downloads", cls="reference internal", href="#downloads"))
            with div(h1("Summary"), cls="section", id="summary"):
                for key,fig in figPlots.items():
                    with div(h2(f"{key}"), cls="section", id=f"{key}"):
                        raw(fig.to_html(full_html=False, include_plotlyjs=PLOTLY_SOURCE))
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
