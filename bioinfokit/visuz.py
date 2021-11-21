"""
visuz module implements visualization functions related to Bioinformatics, Statistics and Machine learning:
Gene expression data visualization
Molecular marker data visualization
Statistical and  Machine learning visualization
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib_venn import venn3, venn2
from random import sample
from functools import reduce
import sys
from matplotlib.colors import ListedColormap

__all__ = ['GeneExpression', 'General', 'gene_exp', 'general', 'marker', 'marker', 'stat', 'cluster']


def venn(vennset=(1,1,1,1,1,1,1), venncolor=('#00909e', '#f67280', '#ff971d'), vennalpha=0.5,
         vennlabel=('A', 'B', 'C')):
    fig = plt.figure()
    if len(vennset) == 7:
        venn3(subsets=vennset, set_labels=vennlabel, set_colors=venncolor, alpha=vennalpha)
        plt.savefig('venn3.png', format='png', bbox_inches='tight', dpi=300)
    elif len(vennset) == 3:
        venn2(subsets=vennset, set_labels=vennlabel, set_colors=venncolor, alpha=vennalpha)
        plt.savefig('venn2.png', format='png', bbox_inches='tight', dpi=300)
    else:
        print("Error: check the set dataset")


class GeneExpression:

    def __init__(self):
        pass

    @staticmethod
    def gene_plot(d, geneid, lfc, lfc_thr, pv_thr, genenames, gfont, pv, gstyle):
        if genenames is not None and genenames == "deg":
            for i in d[geneid].unique():
                if (d.loc[d[geneid] == i, lfc].iloc[0] >= lfc_thr[0] and d.loc[d[geneid] == i, pv].iloc[0] < pv_thr[0]) or \
                        (d.loc[d[geneid] == i, lfc].iloc[0] <= -lfc_thr[1] and d.loc[d[geneid] == i, pv].iloc[0] < pv_thr[1]):
                    if gstyle == 1:
                        plt.text(d.loc[d[geneid] == i, lfc].iloc[0], d.loc[d[geneid] == i, 'logpv_add_axy'].iloc[0], i,
                                      fontsize=gfont)
                    elif gstyle == 2:
                        plt.annotate(i, xy=(d.loc[d[geneid] == i, lfc].iloc[0], d.loc[d[geneid] == i, 'logpv_add_axy'].iloc[0]),
                                     xycoords='data', xytext=(5, -15), textcoords='offset points', size=6,
                                     bbox=dict(boxstyle="round", alpha=0.1),
                                     arrowprops=dict(arrowstyle="wedge,tail_width=0.5", alpha=0.1, relpos=(0, 0)))
                    else:
                        print("Error: invalid gstyle choice")
                        sys.exit(1)
        elif genenames is not None and type(genenames) is tuple:
            for i in d[geneid].unique():
                if i in genenames:
                    if gstyle == 1:
                        plt.text(d.loc[d[geneid] == i, lfc].iloc[0], d.loc[d[geneid] == i, 'logpv_add_axy'].iloc[0], i,
                                      fontsize=gfont)
                    elif gstyle == 2:
                        plt.annotate(i, xy=(d.loc[d[geneid] == i, lfc].iloc[0], d.loc[d[geneid] == i, 'logpv_add_axy'].iloc[0]),
                                     xycoords='data', xytext=(5, -15), textcoords='offset points', size=6,
                                     bbox=dict(boxstyle="round", alpha=0.1),
                                     arrowprops=dict(arrowstyle="wedge,tail_width=0.5", alpha=0.1, relpos=(0, 0)))
                    else:
                        print("Error: invalid gstyle choice")
                        sys.exit(1)
        elif genenames is not None and type(genenames) is dict:
            for i in d[geneid].unique():
                if i in genenames:
                    if gstyle == 1:
                        plt.text(d.loc[d[geneid] == i, lfc].iloc[0], d.loc[d[geneid] == i, 'logpv_add_axy'].iloc[0],
                                      genenames[i], fontsize=gfont)
                    elif gstyle == 2:
                        plt.annotate(genenames[i], xy=(d.loc[d[geneid] == i, lfc].iloc[0], d.loc[d[geneid] == i, 'logpv_add_axy'].iloc[0]),
                                     xycoords='data', xytext=(5, -15), textcoords='offset points', size=6,
                                     bbox=dict(boxstyle="round", alpha=0.1),
                                     arrowprops=dict(arrowstyle="wedge,tail_width=0.5", alpha=0.1, relpos=(0, 0)))
                    else:
                        print("Error: invalid gstyle choice")
                        sys.exit(1)

    @staticmethod
    def geneplot_ma(df, geneid, lfc, lfc_thr, genenames, gfont, gstyle):
        if genenames is not None and genenames == "deg":
            for i in df[geneid].unique():
                if df.loc[df[geneid] == i, lfc].iloc[0] >= lfc_thr[0] or \
                        df.loc[df[geneid] == i, lfc].iloc[0] <= -lfc_thr[1]:
                    if gstyle == 1:
                        plt.text(df.loc[df[geneid] == i, 'A_add_axy'].iloc[0], df.loc[df[geneid] == i, lfc].iloc[0], i,
                                 fontsize=gfont)
                    elif gstyle == 2:
                        plt.annotate(i, xy=(df.loc[df[geneid] == i, 'A_add_axy'].iloc[0],
                                            df.loc[df[geneid] == i, lfc].iloc[0]),
                                     xycoords='data', xytext=(5, -15), textcoords='offset points', size=6,
                                     bbox=dict(boxstyle="round", alpha=0.1),
                                     arrowprops=dict(arrowstyle="wedge,tail_width=0.5", alpha=0.1, relpos=(0, 0)))
                    else:
                        print("Error: invalid gstyle choice")
                        sys.exit(1)
        elif genenames is not None and type(genenames) is tuple:
            for i in df[geneid].unique():
                if i in genenames:
                    if gstyle == 1:
                        plt.text(df.loc[df[geneid] == i, 'A_add_axy'].iloc[0], df.loc[df[geneid] == i, lfc].iloc[0], i,
                                 fontsize=gfont)
                    elif gstyle == 2:
                        plt.annotate(i, xy=(df.loc[df[geneid] == i, 'A_add_axy'].iloc[0],
                                            df.loc[df[geneid] == i, lfc].iloc[0]),
                                     xycoords='data', xytext=(5, -15), textcoords='offset points', size=6,
                                     bbox=dict(boxstyle="round", alpha=0.1),
                                     arrowprops=dict(arrowstyle="wedge,tail_width=0.5", alpha=0.1, relpos=(0, 0)))
                    else:
                        print("Error: invalid gstyle choice")
                        sys.exit(1)
        elif genenames is not None and type(genenames) is dict:
            for i in df[geneid].unique():
                if i in genenames:
                    if gstyle == 1:
                        plt.text(df.loc[df[geneid] == i, 'A_add_axy'].iloc[0], df.loc[df[geneid] == i, lfc].iloc[0],
                                 genenames[i], fontsize=gfont)
                    elif gstyle == 2:
                        plt.annotate(genenames[i], xy=(df.loc[df[geneid] == i, 'A_add_axy'].iloc[0],
                                                       df.loc[df[geneid] == i, lfc].iloc[0]),
                                     xycoords='data', xytext=(5, -15), textcoords='offset points', size=6,
                                     bbox=dict(boxstyle="round", alpha=0.1),
                                     arrowprops=dict(arrowstyle="wedge,tail_width=0.5", alpha=0.1, relpos=(0, 0)))
                    else:
                        print("Error: invalid gstyle choice")
                        sys.exit(1)

    def volcano(df="dataframe", lfc=None, pv=None, lfc_thr=(1, 1), pv_thr=(0.05, 0.05), color=("green", "grey", "red"),
                valpha=1, geneid=None, genenames=None, gfont=8, dim=(5, 5), r=300, ar=90, dotsize=8, markerdot="o",
                sign_line=False, gstyle=1, show=False, figtype='png', axtickfontsize=9,
                axtickfontname="Arial", axlabelfontsize=9, axlabelfontname="Arial", axxlabel=None,
                axylabel=None, xlm=None, ylm=None, plotlegend=False, legendpos='best',
                figname='volcano', legendanchor=None,
                legendlabels=['significant up', 'not significant', 'significant down'], theme=None):
        _x = r'$ log_{2}(Fold Change)$'
        _y = r'$ -log_{10}(P-value)$'
        color = color
        # check if dataframe contains any non-numeric character
        assert general.check_for_nonnumeric(df[lfc]) == 0, 'dataframe contains non-numeric values in lfc column'
        assert general.check_for_nonnumeric(df[pv]) == 0, 'dataframe contains non-numeric values in pv column'
        # this is important to check if color or logpv exists and drop them as if you run multiple times same command
        # it may update old instance of df
        df = df.drop(['color_add_axy', 'logpv_add_axy'], axis=1, errors='ignore')
        assert len(set(color)) == 3, 'unique color must be size of 3'
        df.loc[(df[lfc] >= lfc_thr[0]) & (df[pv] < pv_thr[0]), 'color_add_axy'] = color[0]  # upregulated
        df.loc[(df[lfc] <= -lfc_thr[1]) & (df[pv] < pv_thr[1]), 'color_add_axy'] = color[2]  # downregulated
        df['color_add_axy'].fillna(color[1], inplace=True)  # intermediate
        df['logpv_add_axy'] = -(np.log10(df[pv]))
        # plot
        assign_values = {col: i for i, col in enumerate(color)}
        color_result_num = [assign_values[i] for i in df['color_add_axy']]
        assert len(set(color_result_num)) == 3, \
            'either significant or non-significant genes are missing; try to change lfc_thr or pv_thr to include ' \
            'both significant and non-significant genes'
        if theme == 'dark':
            general.dark_bg()
        plt.subplots(figsize=dim)
        if plotlegend:
            s = plt.scatter(df[lfc], df['logpv_add_axy'], c=color_result_num, cmap=ListedColormap(color), alpha=valpha,
                            s=dotsize, marker=markerdot)
            assert len(legendlabels) == 3, 'legendlabels must be size of 3'
            plt.legend(handles=s.legend_elements()[0], labels=legendlabels, loc=legendpos, bbox_to_anchor=legendanchor)
        else:
            plt.scatter(df[lfc], df['logpv_add_axy'], c=color_result_num, cmap=ListedColormap(color), alpha=valpha,
                        s=dotsize, marker=markerdot)
        if sign_line:
            plt.axhline(y=-np.log10(pv_thr[0]), linestyle='--', color='#7d7d7d', linewidth=1)
            plt.axvline(x=lfc_thr[0], linestyle='--', color='#7d7d7d', linewidth=1)
            plt.axvline(x=-lfc_thr[1], linestyle='--', color='#7d7d7d', linewidth=1)
        GeneExpression.gene_plot(df, geneid, lfc, lfc_thr, pv_thr, genenames, gfont, pv, gstyle)

        if axxlabel:
            _x = axxlabel
        if axylabel:
            _y = axylabel
        general.axis_labels(_x, _y, axlabelfontsize, axlabelfontname)
        general.axis_ticks(xlm, ylm, axtickfontsize, axtickfontname, ar)
        general.get_figure(show, r, figtype, figname, theme)

    def involcano(df="dataframe", lfc="logFC", pv="p_values", lfc_thr=(1, 1), pv_thr=(0.05, 0.05), color=("green", "grey", "red"),
                  valpha=1, geneid=None, genenames=None, gfont=8, dim=(5, 5), r=300, ar=90, dotsize=8, markerdot="o",
                sign_line=False, gstyle=1, show=False, figtype='png', axtickfontsize=9,
               axtickfontname="Arial", axlabelfontsize=9, axlabelfontname="Arial", axxlabel=None,
                axylabel=None, xlm=None, ylm=None, plotlegend=False, legendpos='best',
                figname='involcano', legendanchor=None, legendlabels=['significant up', 'not significant', 'significant down'],
                  theme=None):
        _x = r'$ log_{2}(Fold Change)$'
        _y = r'$ -log_{10}(P-value)$'
        color = color
        assert general.check_for_nonnumeric(df[lfc]) == 0, 'dataframe contains non-numeric values in lfc column'
        assert general.check_for_nonnumeric(df[pv]) == 0, 'dataframe contains non-numeric values in pv column'
        # this is important to check if color or logpv exists and drop them as if you run multiple times same command
        # it may update old instance of df
        df = df.drop(['color_add_axy', 'logpv_add_axy'], axis=1, errors='ignore')
        assert len(set(color)) == 3, 'unique color must be size of 3'
        df.loc[(df[lfc] >= lfc_thr[0]) & (df[pv] < pv_thr[0]), 'color_add_axy'] = color[0]  # upregulated
        df.loc[(df[lfc] <= -lfc_thr[1]) & (df[pv] < pv_thr[1]), 'color_add_axy'] = color[2]  # downregulated
        df['color_add_axy'].fillna(color[1], inplace=True)  # intermediate
        df['logpv_add_axy'] = -(np.log10(df[pv]))

        # plot
        assign_values = {col: i for i, col in enumerate(color)}
        color_result_num = [assign_values[i] for i in df['color_add_axy']]
        assert len(set(color_result_num)) == 3, 'either significant or non-significant genes are missing; try to change lfc_thr or ' \
                                           'pv_thr to include  both significant and non-significant genes'
        if theme == 'dark':
            general.dark_bg()
        plt.subplots(figsize=dim)
        if plotlegend:
            s = plt.scatter(df[lfc], df['logpv_add_axy'], c=color_result_num, cmap=ListedColormap(color), alpha=valpha,
                    s=dotsize, marker=markerdot)
            assert len(legendlabels) == 3, 'legendlabels must be size of 3'
            plt.legend(handles=s.legend_elements()[0], labels=legendlabels, loc=legendpos,
                       bbox_to_anchor=legendanchor)
        else:
            plt.scatter(df[lfc], df['logpv_add_axy'], c=color_result_num, cmap=ListedColormap(color), alpha=valpha,
                        s=dotsize, marker=markerdot)
        GeneExpression.gene_plot(df, geneid, lfc, lfc_thr, pv_thr, genenames, gfont, pv, gstyle)
        plt.gca().invert_yaxis()
        if axxlabel:
            _x = axxlabel
        if axylabel:
            _y = axylabel
        general.axis_labels(_x, _y, axlabelfontsize, axlabelfontname)
        if xlm:
            print('Error: xlm not compatible with involcano')
            sys.exit(1)
        if ylm:
            print('Error: ylm not compatible with involcano')
            sys.exit(1)
        general.axis_ticks(xlm, ylm, axtickfontsize, axtickfontname, ar)
        general.get_figure(show, r, figtype, figname, theme)

    @staticmethod
    def ma(df="dataframe", lfc=None, ct_count=None, st_count=None, basemean=None, pv=None, lfc_thr=(1, 1), pv_thr=0.05,
           valpha=1, dotsize=8,markerdot="o", dim=(6, 5), r=300, show=False, color=("green", "grey", "red"), ar=0,
           figtype='png',axtickfontsize=9, axtickfontname="Arial", axlabelfontsize=9, axlabelfontname="Arial",
           axxlabel=None, axylabel=None, xlm=None, ylm=None, fclines=False, fclinescolor='#2660a4', legendpos='best',
           figname='ma', legendanchor=None, legendlabels=['significant up', 'not significant', 'significant down'],
           plotlegend=False, theme=None, geneid=None, genenames=None, gfont=8, gstyle=1, title=None):
        _x, _y = 'A', 'M'
        assert General.check_for_nonnumeric(df[lfc]) == 0, 'dataframe contains non-numeric values in lfc column'
        if ct_count and st_count:
            assert General.check_for_nonnumeric(df[ct_count]) == 0, \
                'dataframe contains non-numeric values in ct_count column'
            assert General.check_for_nonnumeric(
                df[st_count]) == 0, 'dataframe contains non-numeric values in ct_count column'
        if basemean:
            assert General.check_for_nonnumeric(df[basemean]) == 0, \
                'dataframe contains non-numeric values in basemean column'

        # this is important to check if color or A exists and drop them as if you run multiple times same command
        # it may update old instance of df
        df = df.drop(['color_add_axy', 'A_add_axy'], axis=1, errors='ignore')
        assert len(set(color)) == 3, 'unique color must be size of 3'
        df.loc[(df[lfc] >= lfc_thr[0]) & (df[pv] < pv_thr), 'color_add_axy'] = color[0]  # upregulated
        df.loc[(df[lfc] <= -lfc_thr[1]) & (df[pv] < pv_thr), 'color_add_axy'] = color[2]  # downregulated
        df['color_add_axy'].fillna(color[1], inplace=True)  # intermediate
        if basemean:
            # basemean (mean of normalized counts from DESeq2 results)
            df['A_add_axy'] = df[basemean]
        else:
            df['A_add_axy'] = (np.log2(df[ct_count]) + np.log2(df[st_count])) / 2
        # plot
        assign_values = {col: i for i, col in enumerate(color)}
        color_result_num = [assign_values[i] for i in df['color_add_axy']]
        assert len(
            set(color_result_num)) == 3, 'either significant or non-significant genes are missing; try to change lfc_thr' \
                                         ' to include both significant and non-significant genes'
        if theme:
            General.style_bg(theme)
        plt.subplots(figsize=dim)
        if plotlegend:
            s = plt.scatter(df['A_add_axy'], df[lfc], c=color_result_num, cmap=ListedColormap(color),
                        alpha=valpha, s=dotsize, marker=markerdot)
            assert len(legendlabels) == 3, 'legendlabels must be size of 3'
            plt.legend(handles=s.legend_elements()[0], labels=legendlabels, loc=legendpos,
                           bbox_to_anchor=legendanchor)
        else:
            plt.scatter(df['A_add_axy'], df[lfc], c=color_result_num, cmap=ListedColormap(color),
                        alpha=valpha, s=dotsize, marker=markerdot)
        # draw a central line at M=0
        plt.axhline(y=0, color='#7d7d7d', linestyle='--')
        # draw lfc threshold lines
        if fclines:
            plt.axhline(y=lfc_thr[0], color=fclinescolor, linestyle='--')
            plt.axhline(y=-lfc_thr[1], color=fclinescolor, linestyle='--')
        if axxlabel:
            _x = axxlabel
        if axylabel:
            _y = axylabel
        GeneExpression.geneplot_ma(df, geneid, lfc, lfc_thr, genenames, gfont, gstyle)
        General.axis_labels(_x, _y, axlabelfontsize, axlabelfontname)
        General.axis_ticks(xlm, ylm, axtickfontsize, axtickfontname, ar)
        General.get_figure(show, r, figtype, figname, theme, title)

    @staticmethod
    def hmap(df="dataframe", cmap="seismic", scale=True, dim=(4, 6), rowclus=True, colclus=True, zscore=None, xlabel=True,
             ylabel=True, tickfont=(10, 10), r=300, show=False, figtype='png', figname='heatmap', theme=None):
        # df = df.set_index(d.columns[0])
        # plot heatmap without cluster
        # more cmap: https://matplotlib.org/3.1.0/tutorials/colors/colormaps.html
        if theme == 'dark':
            general.dark_bg()
        fig, hm = plt.subplots(figsize=dim)
        if rowclus and colclus:
            hm = sns.clustermap(df, cmap=cmap, cbar=scale, z_score=zscore, xticklabels=xlabel, yticklabels=ylabel,
                                figsize=dim)
            hm.ax_heatmap.set_xticklabels(hm.ax_heatmap.get_xmajorticklabels(), fontsize=tickfont[0])
            hm.ax_heatmap.set_yticklabels(hm.ax_heatmap.get_ymajorticklabels(), fontsize=tickfont[1])
            general.get_figure(show, r, figtype, figname, theme)
        elif rowclus and colclus is False:
            hm = sns.clustermap(df, cmap=cmap, cbar=scale, z_score=zscore, xticklabels=xlabel, yticklabels=ylabel,
                                figsize=dim, row_cluster=True, col_cluster=False)
            hm.ax_heatmap.set_xticklabels(hm.ax_heatmap.get_xmajorticklabels(), fontsize=tickfont[0])
            hm.ax_heatmap.set_yticklabels(hm.ax_heatmap.get_ymajorticklabels(), fontsize=tickfont[1])
            general.get_figure(show, r, figtype, figname, theme)
        elif colclus and rowclus is False:
            hm = sns.clustermap(df, cmap=cmap, cbar=scale, z_score=zscore, xticklabels=xlabel, yticklabels=ylabel,
                                figsize=dim, row_cluster=False, col_cluster=True)
            hm.ax_heatmap.set_xticklabels(hm.ax_heatmap.get_xmajorticklabels(), fontsize=tickfont[0])
            hm.ax_heatmap.set_yticklabels(hm.ax_heatmap.get_ymajorticklabels(), fontsize=tickfont[1])
            general.get_figure(show, r, figtype, figname, theme)
        else:
            hm = sns.heatmap(df, cmap=cmap, cbar=scale, xticklabels=xlabel, yticklabels=ylabel)
            plt.xticks(fontsize=tickfont[0])
            plt.yticks(fontsize=tickfont[1])
            general.get_figure(show, r, figtype, figname, theme)


class gene_exp:

    def __init__(self):
        pass

    def geneplot(d, geneid, lfc, lfc_thr, pv_thr, genenames, gfont, pv, gstyle):
        if genenames is not None and genenames == "deg":
            for i in d[geneid].unique():
                if (d.loc[d[geneid] == i, lfc].iloc[0] >= lfc_thr[0] and d.loc[d[geneid] == i, pv].iloc[0] < pv_thr[0]) or \
                        (d.loc[d[geneid] == i, lfc].iloc[0] <= -lfc_thr[1] and d.loc[d[geneid] == i, pv].iloc[0] < pv_thr[1]):
                    if gstyle==1:
                        plt.text(d.loc[d[geneid] == i, lfc].iloc[0], d.loc[d[geneid] == i, 'logpv_add_axy'].iloc[0], i,
                                      fontsize=gfont)
                    elif gstyle==2:
                        plt.annotate(i, xy=(d.loc[d[geneid] == i, lfc].iloc[0], d.loc[d[geneid] == i, 'logpv_add_axy'].iloc[0]),
                                     xycoords='data', xytext=(5, -15), textcoords='offset points', size=6,
                                     bbox=dict(boxstyle="round", alpha=0.1),
                                     arrowprops=dict(arrowstyle="wedge,tail_width=0.5", alpha=0.1, relpos=(0, 0)))
                    else:
                        print("Error: invalid gstyle choice")
                        sys.exit(1)
        elif genenames is not None and type(genenames) is tuple:
            for i in d[geneid].unique():
                if i in genenames:
                    if gstyle==1:
                        plt.text(d.loc[d[geneid] == i, lfc].iloc[0], d.loc[d[geneid] == i, 'logpv_add_axy'].iloc[0], i,
                                      fontsize=gfont)
                    elif gstyle==2:
                        plt.annotate(i, xy=(d.loc[d[geneid] == i, lfc].iloc[0], d.loc[d[geneid] == i, 'logpv_add_axy'].iloc[0]),
                                     xycoords='data', xytext=(5, -15), textcoords='offset points', size=6,
                                     bbox=dict(boxstyle="round", alpha=0.1),
                                     arrowprops=dict(arrowstyle="wedge,tail_width=0.5", alpha=0.1, relpos=(0, 0)))
                    else:
                        print("Error: invalid gstyle choice")
                        sys.exit(1)
        elif genenames is not None and type(genenames) is dict:
            for i in d[geneid].unique():
                if i in genenames:
                    if gstyle==1:
                        plt.text(d.loc[d[geneid] == i, lfc].iloc[0], d.loc[d[geneid] == i, 'logpv_add_axy'].iloc[0],
                                      genenames[i], fontsize=gfont)
                    elif gstyle == 2:
                        plt.annotate(genenames[i], xy=(d.loc[d[geneid] == i, lfc].iloc[0], d.loc[d[geneid] == i, 'logpv_add_axy'].iloc[0]),
                                     xycoords='data', xytext=(5, -15), textcoords='offset points', size=6,
                                     bbox=dict(boxstyle="round", alpha=0.1),
                                     arrowprops=dict(arrowstyle="wedge,tail_width=0.5", alpha=0.1, relpos=(0, 0)))
                    else:
                        print("Error: invalid gstyle choice")
                        sys.exit(1)



    def hmap(df="dataframe", cmap="seismic", scale=True, dim=(4, 6), rowclus=True, colclus=True, zscore=None, xlabel=True,
             ylabel=True, tickfont=(10, 10), r=300, show=False, figtype='png', figname='heatmap', theme=None):
        # df = df.set_index(d.columns[0])
        # plot heatmap without cluster
        # more cmap: https://matplotlib.org/3.1.0/tutorials/colors/colormaps.html
        if theme == 'dark':
            general.dark_bg()
        fig, hm = plt.subplots(figsize=dim)
        if rowclus and colclus:
            hm = sns.clustermap(df, cmap=cmap, cbar=scale, z_score=zscore, xticklabels=xlabel, yticklabels=ylabel,
                                figsize=dim)
            hm.ax_heatmap.set_xticklabels(hm.ax_heatmap.get_xmajorticklabels(), fontsize=tickfont[0])
            hm.ax_heatmap.set_yticklabels(hm.ax_heatmap.get_ymajorticklabels(), fontsize=tickfont[1])
            general.get_figure(show, r, figtype, figname, theme)
        elif rowclus and colclus is False:
            hm = sns.clustermap(df, cmap=cmap, cbar=scale, z_score=zscore, xticklabels=xlabel, yticklabels=ylabel,
                                figsize=dim, row_cluster=True, col_cluster=False)
            hm.ax_heatmap.set_xticklabels(hm.ax_heatmap.get_xmajorticklabels(), fontsize=tickfont[0])
            hm.ax_heatmap.set_yticklabels(hm.ax_heatmap.get_ymajorticklabels(), fontsize=tickfont[1])
            general.get_figure(show, r, figtype, figname, theme)
        elif colclus and rowclus is False:
            hm = sns.clustermap(df, cmap=cmap, cbar=scale, z_score=zscore, xticklabels=xlabel, yticklabels=ylabel,
                                figsize=dim, row_cluster=False, col_cluster=True)
            hm.ax_heatmap.set_xticklabels(hm.ax_heatmap.get_xmajorticklabels(), fontsize=tickfont[0])
            hm.ax_heatmap.set_yticklabels(hm.ax_heatmap.get_ymajorticklabels(), fontsize=tickfont[1])
            general.get_figure(show, r, figtype, figname, theme)
        else:
            hm = sns.heatmap(df, cmap=cmap, cbar=scale, xticklabels=xlabel, yticklabels=ylabel)
            plt.xticks(fontsize=tickfont[0])
            plt.yticks(fontsize=tickfont[1])
            general.get_figure(show, r, figtype, figname, theme)


class General:

    rand_colors = ('#a7414a', '#282726', '#6a8a82', '#a37c27', '#563838', '#0584f2', '#f28a30', '#f05837',
                   '#6465a5', '#00743f', '#be9063', '#de8cf0', '#888c46', '#c0334d', '#270101', '#8d2f23',
                   '#ee6c81', '#65734b', '#14325c', '#704307', '#b5b3be', '#f67280', '#ffd082', '#ffd800',
                   '#ad62aa', '#21bf73', '#a0855b', '#5edfff', '#08ffc8', '#ca3e47', '#c9753d', '#6c5ce7')

    def __init__(self):
        pass

    @staticmethod
    def get_figure(show, r, figtype, fig_name, theme, title):
        if title:
            plt.title(title)
        if show:
            plt.show()
        else:
            plt.savefig(fig_name+'.'+figtype, format=figtype, bbox_inches='tight', dpi=r)
        if theme:
            plt.style.use('default')
        plt.clf()
        plt.close()

    @staticmethod
    def axis_labels(x, y, axlabelfontsize=None, axlabelfontname=None):
        plt.xlabel(x, fontsize=axlabelfontsize, fontname=axlabelfontname)
        plt.ylabel(y, fontsize=axlabelfontsize, fontname=axlabelfontname)

    @staticmethod
    def axis_ticks(xlm=None, ylm=None, axtickfontsize=None, axtickfontname=None, ar=None):
        if xlm:
            plt.xlim(left=xlm[0], right=xlm[1])
            plt.xticks(np.arange(xlm[0], xlm[1], xlm[2]), fontsize=axtickfontsize, rotation=ar, fontname=axtickfontname)
        else:
            plt.xticks(fontsize=axtickfontsize, rotation=ar, fontname=axtickfontname)

        if ylm:
            plt.ylim(bottom=ylm[0], top=ylm[1])
            plt.yticks(np.arange(ylm[0], ylm[1], ylm[2]), fontsize=axtickfontsize, rotation=ar, fontname=axtickfontname)
        else:
            plt.yticks(fontsize=axtickfontsize, rotation=ar, fontname=axtickfontname)

    @staticmethod
    def depr_mes(func_name):
        print("This function is deprecated. Please use", func_name)
        print("Read docs at https://reneshbedre.github.io/blog/howtoinstall.html")

    @staticmethod
    def check_for_nonnumeric(pd_series=None):
        if pd.to_numeric(pd_series, errors='coerce').isna().sum() == 0:
            return 0
        else:
            return 1

    @staticmethod
    def pvalue_symbol(pv=None, symbol=None):
        if 0.05 >= pv > 0.01:
            return symbol
        elif 0.01 >= pv > 0.001:
            return 2 * symbol
        elif pv <= 0.001:
            return 3 * symbol
        else:
            return None

    @staticmethod
    def get_file_from_gd(url=None):
        get_path = 'https://drive.google.com/uc?export=download&id=' + url.split('/')[-2]
        return pd.read_csv(get_path, comment='#')

    @staticmethod
    def style_bg(theme=None):
        plt.style.use(theme)


class general:
    def __init__(self):
        pass

    rand_colors = ('#a7414a', '#282726', '#6a8a82', '#a37c27', '#563838', '#0584f2', '#f28a30', '#f05837',
                   '#6465a5', '#00743f', '#be9063', '#de8cf0', '#888c46', '#c0334d', '#270101', '#8d2f23',
                   '#ee6c81', '#65734b', '#14325c', '#704307', '#b5b3be', '#f67280', '#ffd082', '#ffd800',
                   '#ad62aa', '#21bf73', '#a0855b', '#5edfff', '#08ffc8', '#ca3e47', '#c9753d', '#6c5ce7')

    @staticmethod
    def get_figure(show, r, figtype, fig_name, theme):
        if show:
            plt.show()
        else:
            plt.savefig(fig_name+'.'+figtype, format=figtype, bbox_inches='tight', dpi=r)
        if theme == 'dark':
            plt.style.use('default')
        plt.clf()
        plt.close()


    @staticmethod
    def axis_labels(x, y, axlabelfontsize=None, axlabelfontname=None):
        plt.xlabel(x, fontsize=axlabelfontsize, fontname=axlabelfontname)
        plt.ylabel(y, fontsize=axlabelfontsize, fontname=axlabelfontname)
        # plt.xticks(fontsize=9, fontname="sans-serif")
        # plt.yticks(fontsize=9, fontname="sans-serif")

    @staticmethod
    def axis_ticks(xlm=None, ylm=None, axtickfontsize=None, axtickfontname=None, ar=None):
        if xlm:
            plt.xlim(left=xlm[0], right=xlm[1])
            plt.xticks(np.arange(xlm[0], xlm[1], xlm[2]),  fontsize=axtickfontsize, rotation=ar, fontname=axtickfontname)
        else:
            plt.xticks(fontsize=axtickfontsize, rotation=ar, fontname=axtickfontname)

        if ylm:
            plt.ylim(bottom=ylm[0], top=ylm[1])
            plt.yticks(np.arange(ylm[0], ylm[1], ylm[2]),  fontsize=axtickfontsize, rotation=ar, fontname=axtickfontname)
        else:
            plt.yticks(fontsize=axtickfontsize, rotation=ar, fontname=axtickfontname)

    @staticmethod
    def depr_mes(func_name):
        print("This function is deprecated. Please use", func_name )
        print("Read docs at https://reneshbedre.github.io/blog/howtoinstall.html")

    @staticmethod
    def check_for_nonnumeric(pd_series=None):
        if pd.to_numeric(pd_series, errors='coerce').isna().sum() == 0:
            return 0
        else:
            return 1

    @staticmethod
    def pvalue_symbol(pv=None, symbol=None):
        if 0.05 >= pv > 0.01:
            return symbol
        elif 0.01 >= pv > 0.001:
            return 2 * symbol
        elif pv <= 0.001:
            return 3 * symbol
        else:
            return None

    @staticmethod
    def get_file_from_gd(url=None):
        get_path = 'https://drive.google.com/uc?export=download&id=' + url.split('/')[-2]
        return pd.read_csv(get_path, comment='#')

    @staticmethod
    def dark_bg():
        plt.style.use('dark_background')

class marker:

    def __init__(self):
        pass

    def geneplot_mhat(df, markeridcol, chr, pv, gwasp, markernames, gfont, gstyle, ax):
        if markeridcol is not None:
            if markernames is not None and markernames is True:
                for i in df[markeridcol].unique():
                    if df.loc[df[markeridcol] == i, pv].iloc[0] <= gwasp:
                        if gstyle == 1:
                            plt.text(df.loc[df[markeridcol] == i, 'ind'].iloc[0], df.loc[df[markeridcol] == i, 'tpval'].iloc[0],
                                    str(i), fontsize=gfont)
                        elif gstyle == 2:
                            plt.annotate(i, xy=(df.loc[df[markeridcol] == i, 'ind'].iloc[0], df.loc[df[markeridcol] == i, 'tpval'].iloc[0]),
                                         xycoords='data', xytext=(5, -15), textcoords='offset points', size=6,
                                         bbox=dict(boxstyle="round", alpha=0.2),
                                         arrowprops=dict(arrowstyle="wedge,tail_width=0.5", alpha=0.2, relpos=(0, 0)))
            elif markernames is not None and isinstance(markernames, (tuple, list)):
                for i in df[markeridcol].unique():
                    if i in markernames:
                        if gstyle == 1:
                            plt.text(df.loc[df[markeridcol] == i, 'ind'].iloc[0], df.loc[df[markeridcol] == i, 'tpval'].iloc[0],
                                str(i), fontsize=gfont)
                        elif gstyle == 2:
                            plt.annotate(i, xy=(df.loc[df[markeridcol] == i, 'ind'].iloc[0], df.loc[df[markeridcol] == i, 'tpval'].iloc[0]),
                                         xycoords='data', xytext=(5, -15), textcoords='offset points', size=6,
                                         bbox=dict(boxstyle="round", alpha=0.2),
                                         arrowprops=dict(arrowstyle="wedge,tail_width=0.5", alpha=0.2, relpos=(0, 0)))
            elif markernames is not None and isinstance(markernames, dict):
                for i in df[markeridcol].unique():
                    if i in markernames:
                        if gstyle == 1:
                            plt.text(df.loc[df[markeridcol] == i, 'ind'].iloc[0], df.loc[df[markeridcol] == i, 'tpval'].iloc[0],
                                 markernames[i], fontsize=gfont)
                        elif gstyle == 2:
                            plt.annotate(markernames[i], xy=(
                            df.loc[df[markeridcol] == i, 'ind'].iloc[0], df.loc[df[markeridcol] == i, 'tpval'].iloc[0]),
                                         xycoords='data', xytext=(5, -15), textcoords='offset points', size=6,
                                         bbox=dict(boxstyle="round", alpha=0.2),
                                         arrowprops=dict(arrowstyle="wedge,tail_width=0.5", alpha=0.2, relpos=(0, 0)))
        else:
            raise Exception("provide 'markeridcol' parameter")

    def mhat(df="dataframe", chr=None, pv=None, log_scale=True, color=None, dim=(6,4), r=300, ar=90, gwas_sign_line=False,
             gwasp=5E-08, dotsize=8, markeridcol=None, markernames=None, gfont=8, valpha=1, show=False, figtype='png',
             axxlabel=None, axylabel=None, axlabelfontsize=9, axlabelfontname="Arial", axtickfontsize=9,
             axtickfontname="Arial", ylm=None, gstyle=1, figname='manhattan', theme=None):

        _x, _y = 'Chromosomes', r'$ -log_{10}(P)$'
        rand_colors = ('#a7414a', '#282726', '#6a8a82', '#a37c27', '#563838', '#0584f2', '#f28a30', '#f05837',
                       '#6465a5', '#00743f', '#be9063', '#de8cf0', '#888c46', '#c0334d', '#270101', '#8d2f23',
                       '#ee6c81', '#65734b', '#14325c', '#704307', '#b5b3be', '#f67280', '#ffd082', '#ffd800',
                       '#ad62aa', '#21bf73', '#a0855b', '#5edfff', '#08ffc8', '#ca3e47', '#c9753d', '#6c5ce7',
                       '#a997df', '#513b56', '#590925', '#007fff', '#bf1363', '#f39237', '#0a3200', '#8c271e')
        if log_scale:
            # minus log10 of P-value
            df['tpval'] = -np.log10(df[pv])
        else:
            # for Fst values
            df['tpval'] = df[pv]
        # df = df.sort_values(chr)
        # if the column contains numeric strings
        df = df.loc[pd.to_numeric(df[chr], errors='coerce').sort_values().index]
        # add indices
        df['ind'] = range(len(df))
        df_group = df.groupby(chr)
        if color is not None and len(color) == 2:
            color_1 = int(df[chr].nunique() / 2) * [color[0]]
            color_2 = int(df[chr].nunique() / 2) * [color[1]]
            if df[chr].nunique() % 2 == 0:
                color_list = list(reduce(lambda x, y: x+y, zip(color_1, color_2)))
            elif df[chr].nunique() % 2 == 1:
                color_list = list(reduce(lambda x, y: x+y, zip(color_1, color_2)))
                color_list.append(color[0])
        elif color is not None and len(color) == df[chr].nunique():
            color_list = color
        elif color is None:
            # select colors randomly from the list based in number of chr
            color_list = sample(rand_colors, df[chr].nunique())
        else:
            print("Error: in color argument")
            sys.exit(1)

        xlabels = []
        xticks = []
        if theme == 'dark':
            general.dark_bg()
        fig, ax = plt.subplots(figsize=dim)
        i = 0
        for label, df1 in df.groupby(chr):
            df1.plot(kind='scatter', x='ind', y='tpval', color=color_list[i], s=dotsize, alpha=valpha, ax=ax)
            df1_max_ind = df1['ind'].iloc[-1]
            df1_min_ind = df1['ind'].iloc[0]
            xlabels.append(label)
            xticks.append((df1_max_ind - (df1_max_ind - df1_min_ind) / 2))
            i += 1

        # add GWAS significant line
        if gwas_sign_line is True:
            ax.axhline(y=-np.log10(gwasp), linestyle='--', color='#7d7d7d', linewidth=1)
        if markernames is not None:
            marker.geneplot_mhat(df, markeridcol, chr, pv, gwasp, markernames, gfont, gstyle, ax=ax)
        ax.margins(x=0)
        ax.margins(y=0)
        ax.set_xticks(xticks)
        if log_scale:
            ax.set_ylim([0, max(df['tpval'] + 1)])
        if ylm:
            ylm = np.arange(ylm[0], ylm[1], ylm[2])
        else:
            ylm = np.arange(0, max(df['tpval']+1), 1)
        ax.set_yticks(ylm)
        ax.set_xticklabels(xlabels, rotation=ar)
        # ax.set_yticklabels(ylm, fontsize=axtickfontsize, fontname=axtickfontname, rotation=ar)
        if axxlabel:
            _x = axxlabel
        if axylabel:
            _y = axylabel
        ax.set_xlabel(_x, fontsize=axlabelfontsize, fontname=axlabelfontname)
        ax.set_ylabel(_y, fontsize=axlabelfontsize, fontname=axlabelfontname)
        general.get_figure(show, r, figtype, figname, theme)


class Statis:
    def __init__(self):
        pass

    @staticmethod
    def count_plot(df='dataframe', factor=None, dim=(6, 4)):
        # set axis labels to None
        _x = None
        _y = None
        get_factors = df['disease'].value_counts().index
        xbar = np.arange(len(get_factors))
        get_factors_counts = df['disease'].value_counts()


class stat:
    def __init__(self):
        pass

    def bardot(df="dataframe", dim=(6, 4), bw=0.4, colorbar="#f2aa4cff", colordot=["#101820ff"], hbsize=4, r=300, ar=0,
               dotsize=6, valphabar=1, valphadot=1, markerdot="o", errorbar=True, show=False, ylm=None, axtickfontsize=9,
               axtickfontname="Arial", axlabelfontsize=9, axlabelfontname="Arial", yerrlw=None, yerrcw=None, axxlabel=None,
                axylabel=None, figtype='png'):
        # set axis labels to None
        _x = None
        _y = None
        xbar = np.arange(len(df.columns.to_numpy()))
        color_list_bar = colorbar
        color_list_dot = colordot
        if len(color_list_dot) == 1:
            color_list_dot = colordot*len(df.columns.to_numpy())
        if theme == 'dark':
            general.dark_bg()
        plt.subplots(figsize=dim)
        if errorbar:
            plt.bar(x=xbar, height=df.describe().loc['mean'], yerr=df.sem(), width=bw, color=color_list_bar, capsize=hbsize,
                zorder=0, alpha=valphabar, error_kw={'elinewidth': yerrlw, 'capthick': yerrcw})
        else:
            plt.bar(x=xbar, height=df.describe().loc['mean'], width=bw, color=color_list_bar,
                   capsize=hbsize,
                   zorder=0, alpha=valphabar)

        plt.xticks(xbar, df.columns.to_numpy(), fontsize=axtickfontsize, rotation=ar, fontname=axtickfontname)
        if axxlabel:
            _x = axxlabel
        if axylabel:
            _y = axylabel
        general.axis_labels(_x, _y, axlabelfontsize, axlabelfontname)
        # ylm must be tuple of start, end, interval
        if ylm:
            plt.ylim(bottom=ylm[0], top=ylm[1])
            plt.yticks(np.arange(ylm[0], ylm[1], ylm[2]), fontsize=axtickfontsize, fontname=axtickfontname)
        plt.yticks(fontsize=axtickfontsize, rotation=ar, fontname=axtickfontname)
        # add dots
        for cols in range(len(df.columns.to_numpy())):
            # get markers from here https://matplotlib.org/3.1.1/api/markers_api.html
            plt.scatter(x=np.linspace(xbar[cols]-bw/2, xbar[cols]+bw/2, int(df.describe().loc['count'][cols])),
                       y=df[df.columns[cols]].dropna(), s=dotsize, color=color_list_dot[cols], zorder=1, alpha=valphadot,
                       marker=markerdot)
        general.get_figure(show, r, figtype, 'bardot', theme)

    def regplot(df="dataframe", x=None, y=None, yhat=None, dim=(6, 4), colordot='#4a4e4d', colorline='#fe8a71', r=300,
                ar=0, dotsize=6, valphaline=1, valphadot=1, linewidth=1, markerdot="o", show=False, axtickfontsize=9,
               axtickfontname="Arial", axlabelfontsize=9, axlabelfontname="Arial", ylm=None, xlm=None, axxlabel=None,
                axylabel=None, figtype='png', theme=None):
        if theme == 'dark':
            general.dark_bg()
        fig, ax = plt.subplots(figsize=dim)
        plt.scatter(df[x].to_numpy(), df[y].to_numpy(), color=colordot, s=dotsize, alpha=valphadot, marker=markerdot,
                    label='Observed data')
        plt.plot(df[x].to_numpy(), df[yhat].to_numpy(), color=colorline, linewidth=linewidth, alpha=valphaline,
                 label='Regression line')
        if axxlabel:
            x = axxlabel
        if axylabel:
            y = axylabel
        general.axis_labels(x, y, axlabelfontsize, axlabelfontname)
        general.axis_ticks(xlm, ylm, axtickfontsize, axtickfontname, ar)
        plt.legend(fontsize=9)
        general.get_figure(show, r, figtype, 'reg_plot', theme)

    def reg_resid_plot(df="dataframe", yhat=None, resid=None, stdresid=None, dim=(6, 4), colordot='#4a4e4d',
                       colorline='#2ab7ca', r=300, ar=0, dotsize=6, valphaline=1, valphadot=1, linewidth=1,
                       markerdot="o", show=False, figtype='png', theme=None):
        if theme == 'dark':
            general.dark_bg()
        fig, ax = plt.subplots(figsize=dim)
        if resid is not None:
            plt.scatter(df[yhat], df[resid], color=colordot, s=dotsize, alpha=valphadot, marker=markerdot)
            plt.axhline(y=0, color=colorline, linestyle='--', linewidth=linewidth, alpha=valphaline)
            plt.xlabel("Fitted")
            plt.ylabel("Residuals")
            general.get_figure(show, r, figtype, 'resid_plot', theme)
        else:
            print ("Error: Provide residual data")
        if stdresid is not None:
            plt.scatter(df[yhat], df[stdresid], color=colordot, s=dotsize, alpha=valphadot, marker=markerdot)
            plt.axhline(y=0, color=colorline, linestyle='--', linewidth=linewidth, alpha=valphaline)
            plt.xlabel("Fitted")
            plt.ylabel("Standardized Residuals")
            general.get_figure(show, r, figtype, 'std_resid_plot', theme)
        else:
            print ("Error: Provide standardized residual data")

    def corr_mat(df="dataframe", corm="pearson", cmap="seismic", r=300, show=False, dim=(6, 5), axtickfontname="Arial",
                 axtickfontsize=7, ar=90, figtype='png', theme=None):
        if theme == 'dark':
            general.dark_bg()
        d_corr = df.corr(method=corm)
        plt.subplots(figsize=dim)
        plt.matshow(d_corr, vmin=-1, vmax=1, cmap=cmap)
        plt.colorbar()
        cols = list(df)
        ticks = list(range(0, len(list(df))))
        plt.xticks(ticks, cols, fontsize=axtickfontsize, fontname=axtickfontname, rotation=ar)
        plt.yticks(ticks, cols, fontsize=axtickfontsize, fontname=axtickfontname)
        general.get_figure(show, r, figtype, 'corr_mat', theme)

    # for data with pre-calculated mean and SE
    def multi_bar(df="dataframe", dim=(5, 4), colbar=None, colerrorbar=None, bw=0.4, colorbar=None, xbarcol=None, r=300,
                  show=False, axtickfontname="Arial", axtickfontsize=9, ax_x_ticklabel=None, ar=90, figtype='png',
                  figname='multi_bar', valphabar=1, legendpos='best', errorbar=False, yerrlw=None, yerrcw=None,
                  plotlegend=False, hbsize=4, ylm=None, add_sign_line=False, pv=None,
                  sign_line_opts={'symbol': '*', 'fontsize': 8, 'linewidth':0.8, 'arrowstyle': '-', 'dist_y_pos': 2.5,
                                  'dist_y_neg': 4.2}, add_sign_symbol=False, sign_symbol_opts={'symbol': '*',
                                                                                              'fontsize': 8 },
                  dotplot=False, sub_cat=None,
                  sub_cat_opts={'y_neg_dist': 3.5, 'fontsize': 8}, sub_cat_label_dist=None, theme=None):
        xbar = np.arange(df.shape[0])
        xbar_temp = xbar
        if theme == 'dark':
            general.dark_bg()
        fig, ax = plt.subplots(figsize=dim)
        assert len(colbar) >= 2, "number of bar should be atleast 2"
        assert len(colbar) == len(colorbar), "number of color should be equivalent to number of column bars"
        if colbar is not None and isinstance(colbar, (tuple, list)):
            for i in range(len(colbar)):
                if errorbar:
                    ax.bar(x=xbar_temp, height=df[colbar[i]], yerr=df[colerrorbar[i]], width=bw, color=colorbar[i],
                           alpha=valphabar, capsize=hbsize, label=colbar[i], error_kw={'elinewidth': yerrlw,
                                                                                       'capthick': yerrcw})
                    xbar_temp = xbar_temp+bw
                else:
                    ax.bar(x=xbar_temp, height=df[colbar[i]], width=bw, color=colorbar[i], alpha=valphabar,
                           label=colbar[i])
                    xbar_temp = xbar_temp + bw
        ax.set_xticks(xbar+( (bw*(len(colbar)-1)) / (1+(len(colbar)-1)) ))
        if ax_x_ticklabel:
            x_ticklabel = ax_x_ticklabel
        else:
            x_ticklabel = df[xbarcol]
        ax.set_xticklabels(x_ticklabel, fontsize=axtickfontsize, rotation=ar, fontname=axtickfontname)
        # ylm must be tuple of start, end, interval
        if ylm:
            plt.ylim(bottom=ylm[0], top=ylm[1])
            plt.yticks(np.arange(ylm[0], ylm[1], ylm[2]), fontsize=axtickfontsize, fontname=axtickfontname)
        if plotlegend:
            plt.legend(loc=legendpos)

        if dotplot:
            for cols in range(len(df2['factors'].unique())):
                ax.scatter(x=np.linspace(xbar[cols] - bw / 2, xbar[cols] + bw / 2, int(reps)),
                       y=df2[(df2['factors'] == df2['factors'].unique()[cols]) & (df2['sample'] == 'M')]['value'],
                       s=dotsize, color="#7d0013", zorder=1, alpha=valphadot,
                       marker=markerdot)

        if add_sign_line:
            if len(colbar) == 2:
                for i in xbar:
                    x_pos = xbar[i]
                    x_pos_2 = xbar[i] + bw
                    y_pos = df[colbar[0]].to_numpy()[i] + df[colerrorbar[0]].to_numpy()[i]
                    y_pos_2 = df[colbar[1]].to_numpy()[i] + df[colerrorbar[1]].to_numpy()[i]
                    # only if y axis is positive
                    if y_pos > 0:
                        y_pos += 0.5
                        y_pos_2 += 0.5
                        pv_symb = general.pvalue_symbol(pv[i], sign_line_opts['symbol'])
                        if pv_symb:
                            ax.annotate('', xy=(x_pos, y_pos), xytext=(x_pos_2, y_pos),
                                        arrowprops={'connectionstyle': 'bar, armA=50, armB=50, angle=180, fraction=0 ',
                                                    'arrowstyle': sign_line_opts['arrowstyle'],
                                                    'linewidth': sign_line_opts['linewidth']})
                            ax.annotate(pv_symb, xy=(np.mean([x_pos, x_pos_2]),  max(y_pos, y_pos_2) +
                                                     sign_line_opts['dist_y_pos']),
                                        fontsize=sign_line_opts['fontsize'], ha="center")
                    else:
                        y_pos -= 0.5
                        y_pos_2 -= 0.5
                        pv_symb = general.pvalue_symbol(pv[i], sign_line_opts['symbol'])
                        if pv_symb:
                            ax.annotate('', xy=(x_pos, y_pos), xytext=(x_pos_2, y_pos),
                                        arrowprops={'connectionstyle': 'bar, armA=50, armB=50, angle=180, fraction=-1 ',
                                                    'arrowstyle': sign_line_opts['arrowstyle'],
                                                    'linewidth': sign_line_opts['linewidth']})
                            ax.annotate(pv_symb, xy=(np.mean([x_pos, x_pos_2]), min(y_pos_2, y_pos) -
                                                     sign_line_opts['dist_y_neg']),
                                        fontsize=sign_line_opts['fontsize'], ha="center")
        if add_sign_symbol:
            if len(colbar) == 2:
                for i in xbar:
                    x_pos = xbar[i]
                    x_pos_2 = xbar[i] + bw
                    # max value size factor is essential for rel pos of symbol
                    y_pos = df[colbar[0]].to_numpy()[i] + df[colerrorbar[0]].to_numpy()[i] + \
                            (max(df[colbar[0]].to_numpy()) / 20)
                    y_pos_2 = df[colbar[1]].to_numpy()[i] + df[colerrorbar[1]].to_numpy()[i] + \
                              (max(df[colbar[1]].to_numpy()) / 20)
                    # only if y axis is positive
                    if y_pos > 0:
                            pv_symb_1 = general.pvalue_symbol(pv[i][0], sign_symbol_opts['symbol'])
                            pv_symb_2 = general.pvalue_symbol(pv[i][1], sign_symbol_opts['symbol'])
                            if pv_symb_1:
                                plt.annotate(pv_symb_1, xy=(x_pos, y_pos), fontsize=sign_symbol_opts['fontsize'],
                                                ha="center")
                            if pv_symb_2:
                                plt.annotate(pv_symb_2, xy=(x_pos_2, y_pos_2), fontsize=sign_symbol_opts['fontsize'],
                                             ha="center")
            elif len(colbar) == 3:
                for i in xbar:
                    x_pos = xbar[i]
                    x_pos_2 = xbar[i] + bw
                    x_pos_3 = xbar[i] + (2 * bw)
                    # max value size factor is essential for rel pos of symbol
                    y_pos = df[colbar[0]].to_numpy()[i] + df[colerrorbar[0]].to_numpy()[i] + \
                            (max(df[colbar[0]].to_numpy()) / 20)
                    y_pos_2 = df[colbar[1]].to_numpy()[i] + df[colerrorbar[1]].to_numpy()[i] + \
                              (max(df[colbar[1]].to_numpy()) / 20)
                    y_pos_3 = df[colbar[2]].to_numpy()[i] + df[colerrorbar[2]].to_numpy()[i] + \
                              (max(df[colbar[2]].to_numpy()) / 20)
                    # only if y axis is positive
                    if y_pos > 0:
                            pv_symb_1 = general.pvalue_symbol(pv[i][0], sign_symbol_opts['symbol'])
                            pv_symb_2 = general.pvalue_symbol(pv[i][1], sign_symbol_opts['symbol'])
                            pv_symb_3 = general.pvalue_symbol(pv[i][2], sign_symbol_opts['symbol'])
                            if pv_symb_1:
                                plt.annotate(pv_symb_1, xy=(x_pos, y_pos), fontsize=sign_symbol_opts['fontsize'],
                                                ha="center")
                            if pv_symb_2:
                                plt.annotate(pv_symb_2, xy=(x_pos_2, y_pos_2), fontsize=sign_symbol_opts['fontsize'],
                                             ha="center")
                            if pv_symb_3:
                                plt.annotate(pv_symb_3, xy=(x_pos_3, y_pos_3), fontsize=sign_symbol_opts['fontsize'],
                                             ha="center")
        # update this later for min_value
        min_value = 0
        sub_cat_i = 0
        if sub_cat:
            if isinstance(sub_cat, dict):
                for k in sub_cat:
                    if isinstance(k, tuple) and len(k) == 2:
                        cat_x_pos, cat_y_pos, cat_x_pos_2 = k[0], min_value - \
                                                            (sub_cat_opts[
                                                                 'y_neg_dist'] * size_factor_to_start_line), k[1]
                        plt.annotate('', xy=(cat_x_pos - (bw / 2), cat_y_pos),
                                     xytext=(cat_x_pos_2 + (bw / 2), cat_y_pos),
                                     arrowprops={'arrowstyle': '-', 'linewidth': 0.5}, annotation_clip=False)
                        if sub_cat_label_dist and isinstance(sub_cat_label_dist, list):
                            plt.annotate(sub_cat[k], xy=(np.mean([cat_x_pos, cat_x_pos_2]),
                                                         cat_y_pos - size_factor_to_start_line - sub_cat_label_dist[
                                                             sub_cat_i]),
                                         ha="center", fontsize=sub_cat_opts['fontsize'], annotation_clip=False)
                            sub_cat_i += 1
                        else:
                            plt.annotate(sub_cat[k], xy=(np.mean([cat_x_pos, cat_x_pos_2]),
                                                         cat_y_pos - size_factor_to_start_line),
                                         ha="center", fontsize=sub_cat_opts['fontsize'], annotation_clip=False)
                    else:
                        raise KeyError("Sub category keys must be tuple of size 2")

        general.get_figure(show, r, figtype, figname, theme)

    # with replicates values stacked replicates
    # need to work on this later
    def multi_bar_raw(df="dataframe", dim=(5, 4), samp_col_name=None, bw=0.4, colorbar=None, r=300,
                  show=False, axtickfontname="Arial", axtickfontsize=(9, 9), ax_x_ticklabel=None, ar=(0, 90), figtype='png',
                  figname='multi_bar', valphabar=1, legendpos='best', errorbar=False, yerrlw=None, yerrcw=None,
                  plotlegend=False, hbsize=4, ylm=None, add_sign_line=False, pv=None,
                  sign_line_opts={'symbol': '*', 'fontsize': 9, 'linewidth': 0.8, 'arrowstyle': '-', 'dist_y_pos': 2.5,
                                  'dist_y_neg': 4.2}, add_sign_symbol=False,
                      sign_symbol_opts={'symbol': '*', 'fontsize': 9, 'fontname':'Arial', 'rotation':0},
                  dotplot=False, dotplot_opts={'dotsize': 5, 'color':'#7d0013', 'valpha': 1, 'marker': 'o'},
                  sign_line_pairs=None, group_let_df=None, legendanchor=None, legendcols=1, legendfontsize=8,
                  axylabel=None, axxlabel=None, symb_dist=None, axlabelfontsize=(9, 9), axlabelar=(0, 90), sub_cat=None,
                  sub_cat_opts={'y_neg_dist': 3.5, 'fontsize': 9, 'fontname':'Arial'}, sub_cat_label_dist=None,
                      legendlabelframe=False, div_fact=20, legend_columnspacing=None, add_text=None, theme=None):
        if samp_col_name is None or colorbar is None:
            raise ValueError('Invalid value for samp_col_name or colorbar options')
        if theme == 'dark':
            general.dark_bg()
        fig, ax = plt.subplots(figsize=dim)
        sample_list = df[samp_col_name].unique()
        # assert len(sample_list) >= 2, "number of bar should be atleast 2"
        df_mean = df.groupby(samp_col_name).mean().reset_index().set_index(samp_col_name).T
        df_sem = df.groupby(samp_col_name).sem().reset_index().set_index(samp_col_name).T
        colbar = sample_list
        colerrorbar = sample_list
        xbar = np.arange(df_mean.shape[0])
        xbar_temp = xbar
        xbarcol = df_mean.index
        assert len(colbar) == len(colorbar), "number of color should be equivalent to number of column bars"
        df_melt = pd.melt(df.reset_index(), id_vars=[samp_col_name], value_vars=df_mean.index)
        variable_list = df_melt['variable'].unique()
        min_value = (0, min(df_mean.min()))[min(df_mean.min()) < 0]

        if colbar is not None:
            for i in range(len(colbar)):
                if errorbar:
                    ax.bar(x=xbar_temp, height=df_mean[colbar[i]], yerr=df_sem[colerrorbar[i]], width=bw,
                           color=colorbar[i], alpha=valphabar, capsize=hbsize, label=colbar[i],
                           error_kw={'elinewidth': yerrlw, 'capthick': yerrcw})
                    xbar_temp = xbar_temp + bw
                else:
                    ax.bar(x=xbar_temp, height=df_mean[colbar[i]], width=bw, color=colorbar[i], alpha=valphabar,
                           label=colbar[i])
                    xbar_temp = xbar_temp + bw

        bw_fact = bw / 2
        ax.set_xticks(xbar+((len(df_mean.columns)-1) * bw_fact) )
        # ax.set_xticks(xbar + ((bw * (len(colbar) - 1)) / (1 + (len(colbar) - 1))))
        if ax_x_ticklabel:
            x_ticklabel = ax_x_ticklabel
        else:
            x_ticklabel = df[xbarcol]
        ax.set_xticklabels(x_ticklabel, fontsize=axtickfontsize[0], rotation=ar[0], fontname=axtickfontname)
        if axylabel:
            ax.set_ylabel(axylabel, fontsize=axlabelfontsize[1], rotation=axlabelar[1], fontname=axtickfontname)
        if axxlabel:
            ax.set_xlabel(axxlabel, fontsize=axlabelfontsize[0], rotation=axlabelar[0], fontname=axtickfontname)
        # ylm must be tuple of start, end, interval
        if ylm:
            plt.ylim(bottom=ylm[0], top=ylm[1])
            plt.yticks(np.arange(ylm[0], ylm[1], ylm[2]), fontsize=axtickfontsize[1],
                       fontname=axtickfontname)
        if plotlegend:
            plt.legend(loc=legendpos, bbox_to_anchor=legendanchor, ncol=legendcols, fontsize=legendfontsize,
                       frameon=legendlabelframe, columnspacing=legend_columnspacing)

        if isinstance(add_text, list):
            plt.text(add_text[0], add_text[1], add_text[2], fontsize=9, fontfamily='Arial')

        if dotplot:
            for cols in range(len(variable_list)):
                move_fact = 0
                for cols1 in range(len(sample_list)):
                        ax.scatter(x=np.linspace(xbar[cols] - bw_fact + move_fact, xbar[cols] + bw_fact + move_fact,
                                         int(df.groupby(samp_col_name).count().loc[sample_list[cols1], variable_list[cols]])),
                           y=df_melt[(df_melt['variable'] == df_melt['variable'].unique()[cols]) & (
                                       df_melt[samp_col_name] == sample_list[cols1])]['value'], s=dotplot_opts['dotsize'],
                               color=dotplot_opts['color'], zorder=10, alpha=dotplot_opts['valpha'],
                               marker=dotplot_opts['marker'])
                        move_fact += 2 * bw_fact

        size_factor_to_start_line = max(df_mean.max()) / div_fact
        y_pos_dict = dict()
        y_pos_dict_trt = dict()
        if add_sign_line:
            if len(colbar) == 2:
                for i in xbar:
                    x_pos = xbar[i]
                    x_pos_2 = xbar[i] + bw
                    y_pos = df_mean[colbar[0]].to_numpy()[i] + df_sem[colerrorbar[0]].to_numpy()[i]
                    y_pos_2 = df_mean[colbar[1]].to_numpy()[i] + df_sem[colerrorbar[1]].to_numpy()[i]
                    # only if y axis is positive
                    if y_pos > 0:
                        y_pos += 0.5
                        y_pos_2 += 0.5
                        pv_symb = general.pvalue_symbol(pv[i], sign_line_opts['symbol'])
                        if pv_symb:
                            ax.annotate('', xy=(x_pos, max(y_pos, y_pos_2)), xytext=(x_pos_2, max(y_pos, y_pos_2)),
                                        arrowprops={'connectionstyle': 'bar, armA=50, armB=50, angle=180, fraction=0 ',
                                                    'arrowstyle': sign_line_opts['arrowstyle'],
                                                    'linewidth': sign_line_opts['linewidth']})
                            ax.annotate(pv_symb, xy=(np.mean([x_pos, x_pos_2]), max(y_pos, y_pos_2) +
                                                     sign_line_opts['dist_y_pos']),
                                        fontsize=sign_line_opts['fontsize'], ha="center")
                    else:
                        y_pos -= 0.5
                        y_pos_2 -= 0.5
                        pv_symb = general.pvalue_symbol(pv[i], sign_line_opts['symbol'])
                        if pv_symb:
                            ax.annotate('', xy=(x_pos, y_pos), xytext=(x_pos_2, y_pos),
                                        arrowprops={'connectionstyle': 'bar, armA=50, armB=50, angle=180, fraction=-1 ',
                                                    'arrowstyle': sign_line_opts['arrowstyle'],
                                                    'linewidth': sign_line_opts['linewidth']})
                            ax.annotate(pv_symb, xy=(np.mean([x_pos, x_pos_2]), min(y_pos_2, y_pos) -
                                                     sign_line_opts['dist_y_neg']),
                                        fontsize=sign_line_opts['fontsize'], ha="center")
            elif len(colbar) == 3:
                for i in xbar:
                    x_pos = xbar[i]
                    x_pos_2 = xbar[i] + bw
                    x_pos_3 = xbar[i] + (2 * bw)
                    y_pos = df_mean[colbar[0]].to_numpy()[i] + df_sem[colerrorbar[0]].to_numpy()[i]
                    y_pos_2 = df_mean[colbar[1]].to_numpy()[i] + df_sem[colerrorbar[1]].to_numpy()[i]
                    y_pos_3 = df_mean[colbar[2]].to_numpy()[i] + df_sem[colerrorbar[2]].to_numpy()[i]

                    # only if y axis is positive
                    if y_pos > 0:
                        y_pos += size_factor_to_start_line / 2
                        y_pos_2 += size_factor_to_start_line / 2
                        y_pos_3 += size_factor_to_start_line / 2

                        pv_symb1 = general.pvalue_symbol(pv[i][0], sign_line_opts['symbol'])
                        pv_symb2 = general.pvalue_symbol(pv[i][1], sign_line_opts['symbol'])
                        if pv_symb1:
                            if max(y_pos, y_pos_2) >= y_pos_3:
                                pass
                            ax.annotate('', xy=(x_pos, max(y_pos, y_pos_2)), xytext=(x_pos_2, max(y_pos, y_pos_2)),
                                        arrowprops={'connectionstyle': 'bar, armA=50, armB=50, angle=180, fraction=0 ',
                                                    'arrowstyle': sign_line_opts['arrowstyle'],
                                                    'linewidth': sign_line_opts['linewidth']})
                            ax.annotate(pv_symb1, xy=(np.mean([x_pos, x_pos_2]), max(y_pos, y_pos_2) +
                                                     size_factor_to_start_line),
                                        fontsize=sign_line_opts['fontsize'], ha="center")
                        if pv_symb2:
                            if max(y_pos, y_pos_3) < y_pos_2:
                                y_pos_3 = y_pos_2 + (4 * size_factor_to_start_line)
                            ax.annotate('', xy=(x_pos, max(y_pos, y_pos_3)), xytext=(x_pos_3, max(y_pos, y_pos_3)),
                                        arrowprops={'connectionstyle': 'bar, armA=50, armB=50, angle=180, fraction=0 ',
                                                    'arrowstyle': sign_line_opts['arrowstyle'],
                                                    'linewidth': sign_line_opts['linewidth']})
                            ax.annotate(pv_symb2, xy=(np.mean([x_pos, x_pos_3]), max(y_pos, y_pos_3) +
                                                      size_factor_to_start_line),
                                        fontsize=sign_line_opts['fontsize'], ha="center")
                    else:
                        y_pos -= 0.5
                        y_pos_2 -= 0.5
                        pv_symb = general.pvalue_symbol(pv[i], sign_line_opts['symbol'])
                        if pv_symb:
                            ax.annotate('', xy=(x_pos, y_pos), xytext=(x_pos_2, y_pos),
                                        arrowprops={'connectionstyle': 'bar, armA=50, armB=50, angle=180, fraction=-1 ',
                                                    'arrowstyle': sign_line_opts['arrowstyle'],
                                                    'linewidth': sign_line_opts['linewidth']})
                            ax.annotate(pv_symb, xy=(np.mean([x_pos, x_pos_2]), min(y_pos_2, y_pos) -
                                                     sign_line_opts['dist_y_neg']),
                                        fontsize=sign_line_opts['fontsize'], ha="center")

        if add_sign_symbol:
            if len(colbar) == 2:
                for i in xbar:
                    x_pos = xbar[i]
                    x_pos_2 = xbar[i] + bw
                    if symb_dist:
                        # max value size factor is essential for rel pos of symbol
                        y_pos = df_mean[colbar[0]].to_numpy()[i] + df_sem[colerrorbar[0]].to_numpy()[i] + \
                                (max(df_mean[colbar[0]].to_numpy()) / 20) + symb_dist[i][0]
                        y_pos_2 = df_mean[colbar[1]].to_numpy()[i] + df_sem[colerrorbar[1]].to_numpy()[i] + \
                                (max(df_mean[colbar[1]].to_numpy()) / 20) + symb_dist[i][1]
                    else:
                        y_pos = df_mean[colbar[0]].to_numpy()[i] + df_sem[colerrorbar[0]].to_numpy()[i] + \
                                (max(df_mean[colbar[0]].to_numpy()) / 20)
                        y_pos_2 = df_mean[colbar[1]].to_numpy()[i] + df_sem[colerrorbar[1]].to_numpy()[i] + \
                                  (max(df_mean[colbar[1]].to_numpy()) / 20)

                    '''
                    y_pos = df[colbar[0]].to_numpy()[i] + df[colerrorbar[0]].to_numpy()[i] + \
                            (max(df[colbar[0]].to_numpy()) / 20)
                    y_pos_2 = df[colbar[1]].to_numpy()[i] + df[colerrorbar[1]].to_numpy()[i] + \
                              (max(df[colbar[1]].to_numpy()) / 20)
                    '''
                    # group_let_df need index column
                    if isinstance(group_let_df, pd.DataFrame):
                        # only if y axis is positive
                        if y_pos > 0:
                            if not pd.isnull(group_let_df.loc[colbar[0], xbarcol[i]]):
                                plt.annotate(group_let_df.loc[colbar[0], xbarcol[i]], xy=(x_pos, y_pos),
                                             fontsize=sign_symbol_opts['fontsize'], ha='center',
                                             fontfamily=sign_symbol_opts['fontname'],
                                             rotation=sign_symbol_opts['rotation'])
                        if y_pos_2 > 0:
                            if not pd.isnull(group_let_df.loc[colbar[1], xbarcol[i]]):
                                plt.annotate(group_let_df.loc[colbar[1], xbarcol[i]], xy=(x_pos_2, y_pos_2),
                                             fontsize=sign_symbol_opts['fontsize'], ha='center',
                                             fontfamily=sign_symbol_opts['fontname'],
                                             rotation=sign_symbol_opts['rotation'])
                    # only if y axis is positive
                    # need to verify this
                    elif pv:
                        if y_pos > 0:
                            pv_symb_1 = general.pvalue_symbol(pv[i][0], sign_symbol_opts['symbol'])
                            pv_symb_2 = general.pvalue_symbol(pv[i][1], sign_symbol_opts['symbol'])
                            if pv_symb_1:
                                plt.annotate(pv_symb_1, xy=(x_pos, y_pos), fontsize=sign_symbol_opts['fontsize'],
                                         ha="center", fontfamily=sign_symbol_opts['fontname'],
                                             rotation=sign_symbol_opts['rotation'])
                            if pv_symb_2:
                                plt.annotate(pv_symb_2, xy=(x_pos_2, y_pos_2), fontsize=sign_symbol_opts['fontsize'],
                                         ha="center", fontfamily=sign_symbol_opts['fontname'],
                                             rotation=sign_symbol_opts['rotation'])
                    else:
                        raise Exception('Either group dataframe of p value list is required')

            elif len(colbar) == 3:
                for i in xbar:
                    x_pos = xbar[i]
                    x_pos_2 = xbar[i] + bw
                    x_pos_3 = xbar[i] + (2 * bw)

                    if symb_dist:
                        # max value size factor is essential for rel pos of symbol
                        y_pos = df_mean[colbar[0]].to_numpy()[i] + df_sem[colerrorbar[0]].to_numpy()[i] + \
                                (max(df_mean[colbar[0]].to_numpy()) / 20) + symb_dist[i][0]
                        y_pos_2 = df_mean[colbar[1]].to_numpy()[i] + df_sem[colerrorbar[1]].to_numpy()[i] + \
                                (max(df_mean[colbar[1]].to_numpy()) / 20) + symb_dist[i][1]
                        y_pos_3 = df_mean[colbar[2]].to_numpy()[i] + df_sem[colerrorbar[2]].to_numpy()[i] + \
                                  (max(df_mean[colbar[2]].to_numpy()) / 20) + symb_dist[i][2]
                    else:
                        y_pos = df_mean[colbar[0]].to_numpy()[i] + df_sem[colerrorbar[0]].to_numpy()[i] + \
                                (max(df_mean[colbar[0]].to_numpy()) / 20)
                        y_pos_2 = df_mean[colbar[1]].to_numpy()[i] + df_sem[colerrorbar[1]].to_numpy()[i] + \
                                  (max(df_mean[colbar[1]].to_numpy()) / 20)
                        y_pos_3 = df_mean[colbar[2]].to_numpy()[i] + df_sem[colerrorbar[2]].to_numpy()[i] + \
                                  (max(df_mean[colbar[2]].to_numpy()) / 20)

                    # group_let_df need index column
                    if isinstance(group_let_df, pd.DataFrame):
                        if y_pos > 0:
                            plt.annotate(group_let_df.loc[colbar[0], xbarcol[i]], xy=(x_pos, y_pos),
                                         fontsize=sign_symbol_opts['fontsize'], ha="center",
                                         fontfamily=sign_symbol_opts['fontname'], rotation=sign_symbol_opts['rotation'])
                        if y_pos_2 > 0:
                            plt.annotate(group_let_df.loc[colbar[1], xbarcol[i]], xy=(x_pos_2, y_pos_2),
                                         fontsize=sign_symbol_opts['fontsize'], ha="center",
                                         fontfamily=sign_symbol_opts['fontname'], rotation=sign_symbol_opts['rotation'])
                        if y_pos_3 > 0:
                            plt.annotate(group_let_df.loc[colbar[2], xbarcol[i]], xy=(x_pos_3, y_pos_3),
                                         fontsize=sign_symbol_opts['fontsize'], ha="center",
                                         fontfamily=sign_symbol_opts['fontname'], rotation=sign_symbol_opts['rotation'])

                    if pv:
                        # only if y axis is positive
                        if y_pos > 0:
                            pv_symb_1 = general.pvalue_symbol(pv[i][0], sign_symbol_opts['symbol'])
                            pv_symb_2 = general.pvalue_symbol(pv[i][1], sign_symbol_opts['symbol'])
                            pv_symb_3 = general.pvalue_symbol(pv[i][2], sign_symbol_opts['symbol'])
                            if pv_symb_1:
                                plt.annotate(pv_symb_1, xy=(x_pos, y_pos), fontsize=sign_symbol_opts['fontsize'],
                                             ha="center", fontfamily=sign_symbol_opts['fontname'],
                                             rotation=sign_symbol_opts['rotation'])
                            if pv_symb_2:
                                plt.annotate(pv_symb_2, xy=(x_pos_2, y_pos_2), fontsize=sign_symbol_opts['fontsize'],
                                             ha="center", fontfamily=sign_symbol_opts['fontname'],
                                             rotation=sign_symbol_opts['rotation'])
                            if pv_symb_3:
                                plt.annotate(pv_symb_3, xy=(x_pos_3, y_pos_3), fontsize=sign_symbol_opts['fontsize'],
                                             ha="center", fontfamily=sign_symbol_opts['fontname'],
                                             rotation=sign_symbol_opts['rotation'])
            elif len(colbar) == 4:
                for i in xbar:
                    x_pos = xbar[i]
                    x_pos_2 = xbar[i] + bw
                    x_pos_3 = xbar[i] + (2 * bw)
                    x_pos_4 = xbar[i] + (3 * bw)
                    if symb_dist:
                        # max value size factor is essential for rel pos of symbol
                        y_pos = df_mean[colbar[0]].to_numpy()[i] + df_sem[colerrorbar[0]].to_numpy()[i] + \
                                (max(df_mean[colbar[0]].to_numpy()) / 20) + symb_dist[i][0]
                        y_pos_2 = df_mean[colbar[1]].to_numpy()[i] + df_sem[colerrorbar[1]].to_numpy()[i] + \
                                (max(df_mean[colbar[1]].to_numpy()) / 20) + symb_dist[i][1]
                        y_pos_3 = df_mean[colbar[2]].to_numpy()[i] + df_sem[colerrorbar[2]].to_numpy()[i] + \
                                (max(df_mean[colbar[2]].to_numpy()) / 20) + symb_dist[i][2]
                        y_pos_4 = df_mean[colbar[3]].to_numpy()[i] + df_sem[colerrorbar[3]].to_numpy()[i] + \
                                (max(df_mean[colbar[3]].to_numpy()) / 20) + symb_dist[i][3]
                    else:
                        y_pos = df_mean[colbar[0]].to_numpy()[i] + df_sem[colerrorbar[0]].to_numpy()[i] + \
                                (max(df_mean[colbar[0]].to_numpy()) / 20)
                        y_pos_2 = df_mean[colbar[1]].to_numpy()[i] + df_sem[colerrorbar[1]].to_numpy()[i] + \
                                  (max(df_mean[colbar[1]].to_numpy()) / 20)
                        y_pos_3 = df_mean[colbar[2]].to_numpy()[i] + df_sem[colerrorbar[2]].to_numpy()[i] + \
                                  (max(df_mean[colbar[2]].to_numpy()) / 20)
                        y_pos_4 = df_mean[colbar[3]].to_numpy()[i] + df_sem[colerrorbar[3]].to_numpy()[i] + \
                                  (max(df_mean[colbar[3]].to_numpy()) / 20)

                    # group_let_df need index column
                    if isinstance(group_let_df, pd.DataFrame):
                        # only if y axis is positive
                        if y_pos > 0:
                            plt.annotate(group_let_df.loc[colbar[0], xbarcol[i]], xy=(x_pos, y_pos),
                                         fontsize=sign_symbol_opts['fontsize'], ha="center",
                                         fontfamily=sign_symbol_opts['fontname'], rotation=sign_symbol_opts['rotation'])
                        if y_pos_2 > 0:
                            plt.annotate(group_let_df.loc[colbar[1], xbarcol[i]], xy=(x_pos_2, y_pos_2),
                                         fontsize=sign_symbol_opts['fontsize'], ha="center",
                                         fontfamily=sign_symbol_opts['fontname'], rotation=sign_symbol_opts['rotation'])
                        if y_pos_3 > 0:
                            plt.annotate(group_let_df.loc[colbar[2], xbarcol[i]], xy=(x_pos_3, y_pos_3),
                                         fontsize=sign_symbol_opts['fontsize'], ha="center",
                                         fontfamily=sign_symbol_opts['fontname'], rotation=sign_symbol_opts['rotation'])
                        if y_pos_4 > 0:
                            plt.annotate(group_let_df.loc[colbar[3], xbarcol[i]], xy=(x_pos_4, y_pos_4),
                                         fontsize=sign_symbol_opts['fontsize'], ha="center",
                                         fontfamily=sign_symbol_opts['fontname'], rotation=sign_symbol_opts['rotation'])

                    # need to work on this for 4 bars
                    if pv:
                        pv_symb_1 = general.pvalue_symbol(pv[i][0], sign_symbol_opts['symbol'])
                        pv_symb_2 = general.pvalue_symbol(pv[i][1], sign_symbol_opts['symbol'])
                        pv_symb_3 = general.pvalue_symbol(pv[i][2], sign_symbol_opts['symbol'])
                        pv_symb_4 = general.pvalue_symbol(pv[i][3], sign_symbol_opts['symbol'])
                        if pv_symb_1:
                            plt.annotate(pv_symb_1, xy=(x_pos, y_pos), fontsize=sign_symbol_opts['fontsize'],
                                         ha="center", fontfamily=sign_symbol_opts['fontname'],
                                             rotation=sign_symbol_opts['rotation'])
                        if pv_symb_2:
                            plt.annotate(pv_symb_2, xy=(x_pos_2, y_pos_2), fontsize=sign_symbol_opts['fontsize'],
                                         ha="center", fontfamily=sign_symbol_opts['fontname'],
                                             rotation=sign_symbol_opts['rotation'])
                        if pv_symb_3:
                            plt.annotate(pv_symb_3, xy=(x_pos_3, y_pos_3), fontsize=sign_symbol_opts['fontsize'],
                                         ha="center", fontfamily=sign_symbol_opts['fontname'],
                                             rotation=sign_symbol_opts['rotation'])
                        if pv_symb_4:
                            plt.annotate(pv_symb_4, xy=(x_pos_4, y_pos_4), fontsize=sign_symbol_opts['fontsize'],
                                         ha="center", fontfamily=sign_symbol_opts['fontname'],
                                             rotation=sign_symbol_opts['rotation'])
            elif len(colbar) == 5:
                for i in xbar:
                    x_pos = xbar[i]
                    x_pos_2 = xbar[i] + bw
                    x_pos_3 = xbar[i] + (2 * bw)
                    x_pos_4 = xbar[i] + (3 * bw)
                    x_pos_5 = xbar[i] + (4 * bw)
                    # max value size factor is essential for rel pos of symbol
                    if symb_dist:
                        y_pos = df_mean[colbar[0]].to_numpy()[i] + df_sem[colerrorbar[0]].to_numpy()[i] + \
                                (max(df_mean[colbar[0]].to_numpy()) / 20) + symb_dist[i][0]
                        y_pos_2 = df_mean[colbar[1]].to_numpy()[i] + df_sem[colerrorbar[1]].to_numpy()[i] + \
                                  (max(df_mean[colbar[1]].to_numpy()) / 20) + symb_dist[i][1]
                        y_pos_3 = df_mean[colbar[2]].to_numpy()[i] + df_sem[colerrorbar[2]].to_numpy()[i] + \
                                  (max(df_mean[colbar[2]].to_numpy()) / 20) + symb_dist[i][2]
                        y_pos_4 = df_mean[colbar[3]].to_numpy()[i] + df_sem[colerrorbar[3]].to_numpy()[i] + \
                                  (max(df_mean[colbar[3]].to_numpy()) / 20) + symb_dist[i][3]
                        y_pos_5 = df_mean[colbar[4]].to_numpy()[i] + df_sem[colerrorbar[4]].to_numpy()[i] + \
                                  (max(df_mean[colbar[4]].to_numpy()) / 20) + symb_dist[i][4]
                    else:
                        y_pos = df_mean[colbar[0]].to_numpy()[i] + df_sem[colerrorbar[0]].to_numpy()[i] + \
                                (max(df_mean[colbar[0]].to_numpy()) / 20)
                        y_pos_2 = df_mean[colbar[1]].to_numpy()[i] + df_sem[colerrorbar[1]].to_numpy()[i] + \
                                  (max(df_mean[colbar[1]].to_numpy()) / 20)
                        y_pos_3 = df_mean[colbar[2]].to_numpy()[i] + df_sem[colerrorbar[2]].to_numpy()[i] + \
                                  (max(df_mean[colbar[2]].to_numpy()) / 20)
                        y_pos_4 = df_mean[colbar[3]].to_numpy()[i] + df_sem[colerrorbar[3]].to_numpy()[i] + \
                                  (max(df_mean[colbar[3]].to_numpy()) / 20)
                        y_pos_5 = df_mean[colbar[4]].to_numpy()[i] + df_sem[colerrorbar[4]].to_numpy()[i] + \
                                  (max(df_mean[colbar[4]].to_numpy()) / 20)

                    # group_let_df need index column
                    if isinstance(group_let_df, pd.DataFrame):
                        # only if y axis is positive
                        if y_pos > 0:
                            plt.annotate(group_let_df.loc[colbar[0], xbarcol[i]], xy=(x_pos, y_pos),
                                         fontsize=sign_symbol_opts['fontsize'], ha="center")
                        if y_pos_2 > 0:
                            plt.annotate(group_let_df.loc[colbar[1], xbarcol[i]], xy=(x_pos_2, y_pos_2),
                                         fontsize=sign_symbol_opts['fontsize'], ha="center")
                        if y_pos_3 > 0:
                            plt.annotate(group_let_df.loc[colbar[2], xbarcol[i]], xy=(x_pos_3, y_pos_3),
                                         fontsize=sign_symbol_opts['fontsize'], ha="center")
                        if y_pos_4 > 0:
                            plt.annotate(group_let_df.loc[colbar[3], xbarcol[i]], xy=(x_pos_4, y_pos_4),
                                         fontsize=sign_symbol_opts['fontsize'], ha="center")
                        if y_pos_5 > 0:
                            plt.annotate(group_let_df.loc[colbar[4], xbarcol[i]], xy=(x_pos_5, y_pos_5),
                                         fontsize=sign_symbol_opts['fontsize'], ha="center")

                    # need to work on this for 4 bars
                    if pv:
                        pv_symb_1 = general.pvalue_symbol(pv[i][0], sign_symbol_opts['symbol'])
                        pv_symb_2 = general.pvalue_symbol(pv[i][1], sign_symbol_opts['symbol'])
                        pv_symb_3 = general.pvalue_symbol(pv[i][2], sign_symbol_opts['symbol'])
                        if pv_symb_1:
                            plt.annotate(pv_symb_1, xy=(x_pos, y_pos), fontsize=sign_symbol_opts['fontsize'],
                                         ha="center")
                        if pv_symb_2:
                            plt.annotate(pv_symb_2, xy=(x_pos_2, y_pos_2), fontsize=sign_symbol_opts['fontsize'],
                                         ha="center")
                        if pv_symb_3:
                            plt.annotate(pv_symb_3, xy=(x_pos_3, y_pos_3), fontsize=sign_symbol_opts['fontsize'],
                                         ha="center")
        sub_cat_i = 0
        if sub_cat:
            if isinstance(sub_cat, dict):
                for k in sub_cat:
                    if isinstance(k, tuple) and len(k) == 2:
                        cat_x_pos, cat_y_pos, cat_x_pos_2 = k[0], min_value - \
                                                            (sub_cat_opts[
                                                                 'y_neg_dist'] * size_factor_to_start_line), k[1]
                        plt.annotate('', xy=(cat_x_pos - (bw / 2), cat_y_pos),
                                     xytext=(cat_x_pos_2 + (bw / 2), cat_y_pos),
                                     arrowprops={'arrowstyle': '-', 'linewidth': 0.5}, annotation_clip=False)
                        if sub_cat_label_dist and isinstance(sub_cat_label_dist, list):
                            plt.annotate(sub_cat[k], xy=(np.mean([cat_x_pos, cat_x_pos_2]),
                                                         cat_y_pos - size_factor_to_start_line - sub_cat_label_dist[
                                                             sub_cat_i]),
                                         ha="center", fontsize=sub_cat_opts['fontsize'], annotation_clip=False,
                                         fontfamily=sub_cat_opts['fontname'])
                            sub_cat_i += 1
                        else:
                            plt.annotate(sub_cat[k], xy=(np.mean([cat_x_pos, cat_x_pos_2]),
                                                         cat_y_pos - size_factor_to_start_line),
                                         ha="center", fontsize=sub_cat_opts['fontsize'], annotation_clip=False,
                                         fontfamily=sub_cat_opts['fontname'])
                    else:
                        raise KeyError("Sub category keys must be tuple of size 2")

        general.get_figure(show, r, figtype, figname, theme)

    # for data with replicates
    # deprecate dist_y_pos and dist_y_neg (repalce with  size_factor_to_start_line)
    @staticmethod
    def singlebar(df='dataframe', dim=(6, 4), bw=0.4, colorbar='#f2aa4cff', hbsize=4, r=300, ar=(0, 0), valphabar=1,
                  errorbar=True, show=False, ylm=None, axtickfontsize=9, axtickfontname='Arial', ax_x_ticklabel=None,
                  axlabelfontsize=9, axlabelfontname='Arial', yerrlw=None, yerrcw=None, axxlabel=None, axylabel=None,
                  figtype='png', add_sign_line=False, pv=None,
                  sign_line_opts={'symbol': '*', 'fontsize': 9, 'linewidth': 0.5, 'arrowstyle': '-', 'fontname':'Arial'},
                  sign_line_pvals=False,
                  add_sign_symbol=False, sign_symbol_opts={'symbol': '*', 'fontsize': 9, 'rotation':0, 'fontname':'Arial'},
                  sign_line_pairs=None, sub_cat=None, sub_cat_opts={'y_neg_dist': 3.5, 'fontsize': 9, 'fontname':'Arial'},
                  sub_cat_label_dist=None, symb_dist=None, group_let=None, df_format=None, samp_col_name=None,
                  col_order=False, dotplot=False, dotsize=6, colordot=['#101820ff'], valphadot=1, markerdot='o',
                  sign_line_pairs_dist=None, sign_line_pv_symb_dist=None, div_fact=20, add_text=None,
                  figname='singlebar', connectionstyle='bar, armA=50, armB=50, angle=180, fraction=0',
                  std_errs_vis='both', yerrzorder=8, theme=None):
        plt.rcParams['mathtext.fontset'] = 'custom'
        plt.rcParams['mathtext.default'] = 'regular'
        plt.rcParams['mathtext.it'] = 'Arial:italic'
        plt.rcParams['mathtext.bf'] = 'Arial:italic:bold'

        # set axis labels to None
        _x = None
        _y = None
        if df_format == 'stack':
            # sample_list = df[samp_col_name].unique()
            if samp_col_name is None:
                raise ValueError('sample column name required')
            df_mean = df.groupby(samp_col_name).mean().reset_index().set_index(samp_col_name).T
            df_sem = df.groupby(samp_col_name).sem().reset_index().set_index(samp_col_name).T
            if col_order:
                df_mean = df_mean[df[samp_col_name].unique()]
                df_sem = df_sem[df[samp_col_name].unique()]
            bar_h = df_mean.iloc[0]
            bar_se = df_sem.iloc[0]
            sample_list = df_mean.columns.to_numpy()
            # get minimum from df
            min_value = (0, df_mean.iloc[0].min())[df_mean.iloc[0].min() < 0]
        else:
            bar_h = df.describe().loc['mean']
            bar_se = df.sem()
            bar_counts = df.describe().loc['count']
            sample_list = df.columns.to_numpy()
            min_value = (0, min(df.min()))[min(df.min()) < 0]

        if std_errs_vis == 'upper':
            std_errs_vis = [len(bar_se)*[0], bar_se]
        elif std_errs_vis == 'lower':
            std_errs_vis = [bar_se, len(bar_se)*[0]]
        elif std_errs_vis == 'both':
            std_errs_vis = bar_se
        else:
            raise ValueError('In valid value for the std_errs_vis')

        xbar = np.arange(len(sample_list))
        color_list_bar = colorbar
        if theme == 'dark':
            general.dark_bg()
        plt.subplots(figsize=dim)
        if errorbar:
            plt.bar(x=xbar, height=bar_h, yerr=std_errs_vis, width=bw, color=color_list_bar,
                    capsize=hbsize, alpha=valphabar, zorder=5, error_kw={'elinewidth': yerrlw, 'capthick': yerrcw,
                                                                          'zorder': yerrzorder})
        else:
            plt.bar(x=xbar, height=bar_h, width=bw, color=color_list_bar, capsize=hbsize, alpha=valphabar)

        if ax_x_ticklabel:
            x_ticklabel = ax_x_ticklabel
        else:
            x_ticklabel = sample_list

        plt.xticks(ticks=xbar, labels=x_ticklabel, fontsize=axtickfontsize, rotation=ar[0], fontname=axtickfontname)
        if axxlabel:
            _x = axxlabel
        if axylabel:
            _y = axylabel
        general.axis_labels(_x, _y, axlabelfontsize, axlabelfontname)
        # ylm must be tuple of start, end, interval
        if ylm:
            plt.ylim(bottom=ylm[0], top=ylm[1])
            plt.yticks(np.arange(ylm[0], ylm[1], ylm[2]), fontsize=axtickfontsize, fontname=axtickfontname)
        plt.yticks(fontsize=axtickfontsize, rotation=ar[1], fontname=axtickfontname)

        color_list_dot = colordot
        if len(color_list_dot) == 1:
            color_list_dot = colordot * len(sample_list)
        # checked for unstacked data
        if dotplot:
            for cols in range(len(sample_list)):
                plt.scatter(
                    x=np.linspace(xbar[cols] - bw / 2, xbar[cols] + bw / 2, int(bar_counts[cols])),
                    y=df[df.columns[cols]].dropna(), s=dotsize, color=color_list_dot[cols], zorder=10, alpha=valphadot,
                    marker=markerdot)

        size_factor_to_start_line = max(bar_h) / div_fact
        # for only adjacent bars (not for multiple bars with single control)
        if add_sign_line:
            for i in xbar:
                if i % 2 != 0:
                    continue
                x_pos = xbar[i]
                x_pos_2 = xbar[i+1]
                y_pos = df.describe().loc['mean'].to_numpy()[i] + df.sem().to_numpy()[i]
                y_pos_2 = df.describe().loc['mean'].to_numpy()[i+1] + df.sem().to_numpy()[i+1]
                # only if y axis is positive; in future make a function to call it (2 times used)
                if y_pos > 0:
                    y_pos += size_factor_to_start_line
                    y_pos_2 += size_factor_to_start_line
                    pv_symb = general.pvalue_symbol(pv[int(i/2)], sign_line_opts['symbol'])
                    if pv_symb:
                        plt.annotate('', xy=(x_pos, max(y_pos, y_pos_2)), xytext=(x_pos_2, max(y_pos, y_pos_2)),
                                    arrowprops={'connectionstyle': connectionstyle,
                                                'arrowstyle': sign_line_opts['arrowstyle'],
                                                'linewidth': sign_line_opts['linewidth']})
                        plt.annotate(pv_symb, xy=(np.mean([x_pos, x_pos_2]), max(y_pos, y_pos_2) +
                                                 sign_line_opts['dist_y_pos']),
                                    fontsize=sign_line_opts['fontsize'], ha="center")

        # for only adjacent bars with one control but multiple treatments
        # need to work for sign_line_pairs (update df on line 1276)
        p_index = 0
        y_pos_dict = dict()
        y_pos_dict_trt = dict()
        if sign_line_pairs:
            for i in sign_line_pairs:
                y_pos_adj = 0
                x_pos = xbar[i[0]]
                x_pos_2 = xbar[i[1]]
                y_pos = df.describe().loc['mean'].to_numpy()[i[0]] + df.sem().to_numpy()[i[0]]
                y_pos_2 = df.describe().loc['mean'].to_numpy()[i[1]] + df.sem().to_numpy()[i[1]]
                # only if y axis is positive; in future make a function to call it (2 times used)
                if y_pos > 0:
                    y_pos += size_factor_to_start_line/2
                    y_pos_2 += size_factor_to_start_line/2
                    # check if the mean of y_pos is not lesser than not other treatments which lies between
                    # eg if 0-1 has higher sign bar than the 0-2
                    if i[0] in y_pos_dict_trt:
                        y_pos_adj = 1
                        if y_pos_2 <= y_pos_dict_trt[i[0]][1]:
                            if sign_line_pairs_dist:
                                y_pos_2 += (y_pos_dict_trt[i[0]][1] - y_pos_2) + (3 * size_factor_to_start_line) + \
                                       sign_line_pairs_dist[p_index]
                            else:
                                y_pos_2 += (y_pos_dict_trt[i[0]][1] - y_pos_2) + (3 * size_factor_to_start_line)
                        elif y_pos <= y_pos_dict_trt[i[0]][0]:
                            if sign_line_pairs_dist:
                                y_pos += 3 * size_factor_to_start_line + sign_line_pairs_dist[p_index]
                            else:
                                y_pos += 3 * size_factor_to_start_line
                    # check if difference is not equivalent between two y_pos
                    # if yes add some distance, so that sign bar will not overlap
                    if i[0] in y_pos_dict:
                        y_pos_adj = 1
                        if 0.75 < df.describe().loc['mean'].to_numpy()[i[0]]/df.describe().loc['mean'].to_numpy()[i[1]] < 1.25:
                            if sign_line_pairs_dist:
                                y_pos += 2 * size_factor_to_start_line + sign_line_pairs_dist[p_index]
                            else:
                                y_pos += 2 * size_factor_to_start_line

                    if y_pos_adj == 0 and sign_line_pairs_dist:
                        if y_pos >= y_pos_2:
                            y_pos += sign_line_pairs_dist[p_index]
                        else:
                            y_pos_2 += sign_line_pairs_dist[p_index]

                    # sign_line_pvals passed, used p values instead of symbols
                    if sign_line_pvals:
                        pv_symb = '$\it{p}$'+ str(pv[p_index])
                    else:
                        pv_symb = general.pvalue_symbol(pv[p_index], sign_line_opts['symbol'])
                    y_pos_dict[i[0]] = y_pos
                    y_pos_dict_trt[i[0]] = [y_pos, y_pos_2]
                    if pv_symb:
                        plt.annotate('', xy=(x_pos, max(y_pos, y_pos_2)), xytext=(x_pos_2, max(y_pos, y_pos_2)),
                                     arrowprops={'connectionstyle': connectionstyle,
                                                     'arrowstyle': sign_line_opts['arrowstyle'],
                                                     'linewidth': sign_line_opts['linewidth']})
                        # here size factor size_factor_to_start_line added instead of sign_line_opts['dist_y_pos']
                        # make this change everywhere in future release
                        plt.annotate(pv_symb, xy=(np.mean([x_pos, x_pos_2]), max(y_pos, y_pos_2) +
                                                  size_factor_to_start_line + sign_line_pv_symb_dist[p_index]),
                                     fontsize=sign_line_opts['fontsize'], ha="center")
                    p_index += 1

        if add_sign_symbol:
            for i in xbar:
                x_pos = xbar[i]
                # y_pos = df.describe().loc['mean'].to_numpy()[i] + df.sem().to_numpy()[i] + size_factor_to_start_line

                if symb_dist:
                    y_pos = bar_h.to_numpy()[i] + bar_se.to_numpy()[i] + \
                            size_factor_to_start_line + symb_dist[i]
                else:
                    y_pos = bar_h.to_numpy()[i] + bar_se.to_numpy()[i] + \
                            size_factor_to_start_line

                # group_let list
                if isinstance(group_let, list):
                    if y_pos > 0:
                        plt.annotate(group_let[i], xy=(x_pos, y_pos),
                                     fontsize=sign_symbol_opts['fontsize'], ha="center",
                                     rotation=sign_symbol_opts['rotation'], fontfamily=sign_symbol_opts['fontname'])

                # only if y axis is positive
                if pv:
                    if y_pos > 0:
                        pv_symb = general.pvalue_symbol(pv[i], sign_symbol_opts['symbol'])
                        if pv_symb:
                            plt.annotate(pv_symb, xy=(x_pos, y_pos), fontsize=sign_symbol_opts['fontsize'], ha="center",
                                         rotation=sign_symbol_opts['rotation'], fontfamily=sign_symbol_opts['fontname'])

        sub_cat_i = 0
        if sub_cat:
            if isinstance(sub_cat, dict):
                for k in sub_cat:
                    if isinstance(k, tuple) and len(k) == 2:
                        cat_x_pos, cat_y_pos, cat_x_pos_2 = k[0], min_value - \
                                                            (sub_cat_opts['y_neg_dist']*size_factor_to_start_line), k[1]
                        plt.annotate('', xy=(cat_x_pos-(bw/2), cat_y_pos), xytext=(cat_x_pos_2+(bw/2), cat_y_pos),
                                     arrowprops={'arrowstyle': '-', 'linewidth': 0.5}, annotation_clip=False)
                        if sub_cat_label_dist and isinstance(sub_cat_label_dist, list):
                            plt.annotate(sub_cat[k], xy=(np.mean([cat_x_pos, cat_x_pos_2]),
                                                         cat_y_pos - size_factor_to_start_line - sub_cat_label_dist[sub_cat_i]),
                                         ha="center", fontsize=sub_cat_opts['fontsize'], annotation_clip=False,
                                         fontfamily=sub_cat_opts['fontname'])
                            sub_cat_i += 1
                        else:
                            plt.annotate(sub_cat[k], xy=(np.mean([cat_x_pos, cat_x_pos_2]),
                                                     cat_y_pos-size_factor_to_start_line),
                                     ha="center", fontsize=sub_cat_opts['fontsize'], annotation_clip=False,
                                         fontfamily=sub_cat_opts['fontname'])
                    else:
                        raise KeyError("Sub category keys must be tuple of size 2")

        if isinstance(add_text, list):
            plt.text(add_text[0], add_text[1], add_text[2], fontsize=9, fontfamily='Arial')

        general.get_figure(show, r, figtype, figname, theme)

    @staticmethod
    def normal_bar(df='dataframe', x_col_name=None, y_col_name=None, dim=(6, 4), bw=0.4, colorbar="#f2aa4cff", r=300,
                   ar=(0, 0), valphabar=1, show=False, ylm=None, axtickfontsize=9, axtickfontname='Arial',
                   ax_x_ticklabel=None, axlabelfontsize=9, axlabelfontname='Arial', axxlabel=None, axylabel=None,
                   figtype='png', figname='normal_bar', theme=None):
        # set axis labels to None
        _x = None
        _y = None
        xbar = np.arange(len(df[x_col_name]))
        if theme == 'dark':
            general.dark_bg()
        plt.subplots(figsize=dim)
        plt.bar(x=xbar, height=df[y_col_name], width=bw, color=colorbar, alpha=valphabar)
        if ax_x_ticklabel:
            x_ticklabel = ax_x_ticklabel
        else:
            x_ticklabel = df[x_col_name].to_numpy()
        plt.xticks(ticks=xbar, labels=x_ticklabel, fontsize=axtickfontsize, rotation=ar[0], fontname=axtickfontname)
        if axxlabel:
            _x = axxlabel
        if axylabel:
            _y = axylabel
        general.axis_labels(_x, _y, axlabelfontsize, axlabelfontname)
        general.get_figure(show, r, figtype, figname, theme)

    def boxplot_single_factor(df='dataframe', column_names=None, grid=False, ar=(0, 0), axtickfontsize=9,
                              axtickfontname='Arial', dim=(6, 4), show=False, figtype='png', figname='boxplot', r=300,
                              ylm=None, box_line_style='-', box_line_width=1, box_line_color='b', med_line_style='-',
                              med_line_width=1, med_line_color='g', whisk_line_color='b', cap_color='b',
                              add_sign_symbol=False, symb_dist=None, sign_symbol_opts={'symbol': '*', 'fontsize': 8 },
                              pv=None, notch=False, outliers=True, fill_box_color=True, dotplot=False, dotsize=6,
                              colordot=['#101820ff'], valphadot=1, markerdot='o', theme=None):
        if theme == 'dark':
            general.dark_bg()
        plt.subplots()
        if column_names:
            xbar = column_names
        else:
            xbar = list(df.columns)
            # rot is x axis rotation
        other_args = {'grid': grid, 'rot': ar[0], 'fontsize': axtickfontsize, 'notch':notch, 'showfliers':outliers,
                      'figsize': dim, 'patch_artist': fill_box_color}
        color_args = {'medians': med_line_color, 'boxes': box_line_color, 'whiskers': whisk_line_color,
                      'caps': cap_color}
        medianprops_args = {'linestyle': med_line_style, 'linewidth': med_line_width}
        boxprops_args = {'linestyle': box_line_style, 'linewidth': box_line_width}

        if isinstance(column_names, list):
            df.boxplot(column=column_names, **other_args, boxprops=boxprops_args, medianprops=medianprops_args,
                       color=color_args)
        else:
            df.boxplot(**other_args, boxprops=boxprops_args, color=color_args, medianprops=medianprops_args)

        # ylm must be tuple of start, end, interval
        if ylm:
            plt.ylim(bottom=ylm[0], top=ylm[1])
            plt.yticks(np.arange(ylm[0], ylm[1], ylm[2]), fontsize=axtickfontsize, fontname=axtickfontname)
        plt.yticks(fontsize=axtickfontsize, rotation=ar[1], fontname=axtickfontname)

        color_list_dot = colordot
        if len(color_list_dot) == 1:
            color_list_dot = colordot * len(xbar)
        # checked for unstacked data
        if dotplot:
            for cols in range(len(xbar)):
                plt.scatter(
                    x=np.linspace(xbar[cols] - bw / 2, xbar[cols] + bw / 2, int(bar_counts[cols])),
                    y=df[df.columns[cols]].dropna(), s=dotsize, color=color_list_dot[cols], zorder=10, alpha=valphadot,
                    marker=markerdot)

        size_factor_to_start_line = max(df.max()) / 20
        if add_sign_symbol:
            # p and symb_dist should be dict
            if isinstance(pv, dict):
                for k, v in pv.items():
                    if isinstance(symb_dist, dict):
                        if k not in symb_dist:
                            symb_dist[k] = 0
                        y_pos = df[k].max() + size_factor_to_start_line + symb_dist[k]
                    else:
                        y_pos = df[k].max() + size_factor_to_start_line

                    if y_pos > 0 and v <= 0.05:
                        pv_symb = general.pvalue_symbol(v, sign_symbol_opts['symbol'])
                        if pv_symb:
                            plt.annotate(pv_symb, xy=((xbar.index(k))+1, y_pos),
                                         fontsize=sign_symbol_opts['fontsize'],
                                         ha="center")

        general.get_figure(show, r, figtype, figname, theme)

    @staticmethod
    def roc(fpr=None, tpr=None, c_line_style='-', c_line_color='#f05f21', c_line_width=1, diag_line=True,
            diag_line_style='--', diag_line_width=1, diag_line_color='b', auc=None, shade_auc=False,
            shade_auc_color='#f48d60',
            axxlabel='False Positive Rate (1 - Specificity)', axylabel='True Positive Rate (Sensitivity)', ar=(0, 0),
            axtickfontsize=9, axtickfontname='Arial', axlabelfontsize=9, axlabelfontname='Arial',
            plotlegend=True, legendpos='lower right', legendanchor=None, legendcols=1, legendfontsize=8,
            legendlabelframe=False, legend_columnspacing=None, per_class=False, dim=(6, 5), show=False, figtype='png',
            figname='roc', r=300, ylm=None, theme=None):
        if theme == 'dark':
            general.dark_bg()
        plt.subplots(figsize=dim)
        # plt.margins(x=0)
        if auc:
            plt.plot(fpr, tpr, color=c_line_color, linestyle=c_line_style, linewidth=c_line_width,
                 label='AUC = %0.4f' % auc)
        else:
            plt.plot(fpr, tpr, color=c_line_color, linestyle=c_line_style, linewidth=c_line_width)
        if diag_line:
            plt.plot([0, 1], [0, 1], color=diag_line_color, linestyle=diag_line_style, linewidth=diag_line_width,
                     label='Chance level')
        if per_class:
            plt.plot([0, 0], [0, 1], color='grey', linestyle='-', linewidth=1)
            plt.plot([0, 1], [1, 1], color='grey', linestyle='-', linewidth=1, label='Perfect performance')
        # ylm must be tuple of start, end, interval
        if ylm:
            plt.ylim(bottom=ylm[0], top=ylm[1])
            plt.yticks(np.arange(ylm[0], ylm[1], ylm[2]), fontsize=axtickfontsize, fontname=axtickfontname)
        plt.yticks(fontsize=axtickfontsize, rotation=ar[1], fontname=axtickfontname)
        if axxlabel:
            _x = axxlabel
        if axylabel:
            _y = axylabel
        if shade_auc:
            plt.fill_between(x=fpr, y1=tpr, color=shade_auc_color)
        if plotlegend:
            plt.legend(loc=legendpos, bbox_to_anchor=legendanchor, ncol=legendcols, fontsize=legendfontsize,
                       frameon=legendlabelframe, columnspacing=legend_columnspacing)
        general.axis_labels(_x, _y, axlabelfontsize, axlabelfontname)
        general.get_figure(show, r, figtype, figname, theme)


class cluster:
    def __init__(self):
        pass

    @staticmethod
    def screeplot(obj="pcascree", axlabelfontsize=9, axlabelfontname="Arial", axxlabel=None,
                axylabel=None, figtype='png', r=300, show=False, dim=(6, 4), theme=None):
        if theme == 'dark':
            general.dark_bg()
        y = [x * 100 for x in obj[1]]
        plt.subplots(figsize=dim)
        plt.bar(obj[0], y)
        xlab='PCs'
        ylab='Proportion of variance (%)'
        if axxlabel:
            xlab = axxlabel
        if axylabel:
            ylab = axylabel
        plt.xticks(fontsize=7, rotation=70)
        general.axis_labels(xlab, ylab, axlabelfontsize, axlabelfontname)
        general.get_figure(show, r, figtype, 'screeplot', theme)

    @staticmethod
    def pcaplot(x=None, y=None, z=None, labels=None, var1=None, var2=None, var3=None, axlabelfontsize=9,
                axlabelfontname="Arial", figtype='png', r=300, show=False, plotlabels=True, dim=(6, 4), theme=None):
        if theme == 'dark':
            general.dark_bg()
        if x is not None and y is not None and z is None:
            assert var1 is not None and var2 is not None and labels is not None, "var1 or var2 variable or labels are missing"
            plt.subplots(figsize=dim)
            for i, varnames in enumerate(labels):
                plt.scatter(x[i], y[i])
                if plotlabels:
                    plt.text(x[i], y[i], varnames, fontsize=10)
            general.axis_labels("PC1 ({}%)".format(var1), "PC2 ({}%)".format(var2), axlabelfontsize, axlabelfontname)
            general.get_figure(show, r, figtype, 'pcaplot_2d', theme)
        elif x is not None and y is not None and z is not None:
            assert var1 and var2 and var3 and labels is not None, "var1 or var2 or var3 or labels are missing"
            # for 3d plot
            fig = plt.figure(figsize=dim)
            ax = fig.add_subplot(111, projection='3d')
            for i, varnames in enumerate(labels):
                ax.scatter(x[i], y[i], z[i])
                if plotlabels:
                    ax.text(x[i], y[i], z[i], varnames, fontsize=10)
            ax.set_xlabel("PC1 ({}%)".format(var1), fontsize=axlabelfontsize, fontname=axlabelfontname)
            ax.set_ylabel("PC2 ({}%)".format(var2), fontsize=axlabelfontsize, fontname=axlabelfontname)
            ax.set_zlabel("PC3 ({}%)".format(var3), fontsize=axlabelfontsize, fontname=axlabelfontname)
            general.get_figure(show, r, figtype, 'pcaplot_3d', theme)

    @staticmethod
    # adapted from https://stackoverflow.com/questions/39216897/plot-pca-loadings-and-loading-in-biplot-in-sklearn-like-rs-autoplot
    def biplot(cscore=None, loadings=None, labels=None, var1=None, var2=None, var3=None, axlabelfontsize=9, axlabelfontname="Arial",
               figtype='png', r=300, show=False, markerdot="o", dotsize=6, valphadot=1, colordot='#eba487', arrowcolor='#87ceeb',
               valphaarrow=1, arrowlinestyle='-', arrowlinewidth=0.5, centerlines=True, colorlist=None, legendpos='best',
               datapoints=True, dim=(6, 4), theme=None):
        if theme == 'dark':
            general.dark_bg()
        assert cscore is not None and loadings is not None and labels is not None and var1 is not None and var2 is not None, \
            "cscore or loadings or labels or var1 or var2 are missing"
        if var1 is not None and var2 is not None and var3 is None:
            xscale = 1.0 / (cscore[:, 0].max() - cscore[:, 0].min())
            yscale = 1.0 / (cscore[:, 1].max() - cscore[:, 1].min())
            # zscale = 1.0 / (cscore[:, 2].max() - cscore[:, 2].min())
            # colorlist is an array of classes from dataframe column
            plt.subplots(figsize=dim)
            if datapoints:
                if colorlist is not None:
                    unique_class = set(colorlist)
                    # color_dict = dict()
                    assign_values = {col: i for i, col in enumerate(unique_class)}
                    color_result_num = [assign_values[i] for i in colorlist]
                    if colordot and isinstance(colordot, (tuple, list)):
                        colour_map = ListedColormap(colordot)
                        # for i in range(len(list(unique_class))):
                        #    color_dict[list(unique_class)[i]] = colordot[i]
                        # color_result = [color_dict[i] for i in colorlist]
                        s = plt.scatter(cscore[:, 0] * xscale, cscore[:, 1] * yscale, c=color_result_num, cmap=colour_map,
                                        s=dotsize, alpha=valphadot, marker=markerdot)
                        plt.legend(handles=s.legend_elements()[0], labels=list(unique_class), loc=legendpos)
                    elif colordot and not isinstance(colordot, (tuple, list)):
                        # s = plt.scatter(cscore[:, 0] * xscale, cscore[:, 1] * yscale, color=color_result, s=dotsize,
                        #                alpha=valphadot, marker=markerdot)
                        # plt.legend(handles=s.legend_elements()[0], labels=list(unique_class))
                        s = plt.scatter(cscore[:, 0] * xscale, cscore[:, 1] * yscale, c=color_result_num, s=dotsize,
                                    alpha=valphadot, marker=markerdot)
                        plt.legend(handles=s.legend_elements()[0], labels=list(unique_class), loc=legendpos)
                else:
                    plt.scatter(cscore[:, 0] * xscale, cscore[:, 1] * yscale, color=colordot, s=dotsize,
                                    alpha=valphadot, marker=markerdot)
            if centerlines:
                plt.axhline(y=0, linestyle='--', color='#7d7d7d', linewidth=1)
                plt.axvline(x=0, linestyle='--', color='#7d7d7d', linewidth=1)
            # loadings[0] is the number of the original variables
            # this is important where variables more than number of observations
            for i in range(len(loadings[0])):
                plt.arrow(0, 0, loadings[0][i], loadings[1][i], color=arrowcolor, alpha=valphaarrow, ls=arrowlinestyle,
                          lw=arrowlinewidth)
                plt.text(loadings[0][i], loadings[1][i], labels[i])
                # adjust_text(t)
            # plt.xlim(min(loadings[0]) - 0.1, max(loadings[0]) + 0.1)
            # plt.ylim(min(loadings[1]) - 0.1, max(loadings[1]) + 0.1)
            xlimit_max = np.max([np.max(cscore[:, 0]*xscale), np.max(loadings[0])])
            xlimit_min = np.min([np.min(cscore[:, 0]*xscale), np.min(loadings[0])])
            ylimit_max = np.max([np.max(cscore[:, 1]*yscale), np.max(loadings[1])])
            ylimit_min = np.min([np.min(cscore[:, 1]*yscale), np.min(loadings[1])])
            plt.xlim(xlimit_min-0.2, xlimit_max+0.2)
            plt.ylim(ylimit_min-0.2, ylimit_max+0.2)
            general.axis_labels("PC1 ({}%)".format(var1), "PC2 ({}%)".format(var2), axlabelfontsize, axlabelfontname)
            general.get_figure(show, r, figtype, 'biplot_2d', theme)
        # 3D
        if var1 is not None and var2 is not None and var3 is not None:
            xscale = 1.0 / (cscore[:, 0].max() - cscore[:, 0].min())
            yscale = 1.0 / (cscore[:, 1].max() - cscore[:, 1].min())
            zscale = 1.0 / (cscore[:, 2].max() - cscore[:, 2].min())
            fig = plt.figure(figsize=dim)
            ax = fig.add_subplot(111, projection='3d')
            if datapoints:
                if colorlist is not None:
                    unique_class = set(colorlist)
                    assign_values = {col: i for i, col in enumerate(unique_class)}
                    color_result_num = [assign_values[i] for i in colorlist]
                    if colordot and isinstance(colordot, (tuple, list)):
                        colour_map = ListedColormap(colordot)
                        s = ax.scatter(cscore[:, 0]*xscale, cscore[:, 1]*yscale, cscore[:, 2]*zscale, c=color_result_num,
                                       cmap=colour_map, s=dotsize, alpha=valphadot, marker=markerdot)
                        plt.legend(handles=s.legend_elements()[0], labels=list(unique_class), loc=legendpos)
                    elif colordot and not isinstance(colordot, (tuple, list)):
                        s = ax.scatter(cscore[:, 0]*xscale, cscore[:, 1]*yscale, cscore[:, 2]*zscale, c=color_result_num,
                                        s=dotsize, alpha=valphadot, marker=markerdot)
                        plt.legend(handles=s.legend_elements()[0], labels=list(unique_class), loc=legendpos)
                else:
                    ax.scatter(cscore[:, 0] * xscale, cscore[:, 1] * yscale, cscore[:, 2] * zscale, color=colordot,
                               s=dotsize, alpha=valphadot, marker=markerdot)
            for i in range(len(loadings[0])):
                ax.quiver(0, 0, 0, loadings[0][i], loadings[1][i], loadings[2][i], color=arrowcolor, alpha=valphaarrow,
                          ls=arrowlinestyle, lw=arrowlinewidth)
                ax.text(loadings[0][i], loadings[1][i], loadings[2][i],  labels[i])

            xlimit_max = np.max([np.max(cscore[:, 0] * xscale), np.max(loadings[0])])

            xlimit_min = np.min([np.min(cscore[:, 0] * xscale), np.min(loadings[0])])
            ylimit_max = np.max([np.max(cscore[:, 1] * yscale), np.max(loadings[1])])
            ylimit_min = np.min([np.min(cscore[:, 1] * yscale), np.min(loadings[1])])
            zlimit_max = np.max([np.max(cscore[:, 2] * zscale), np.max(loadings[2])])
            zlimit_min = np.min([np.min(cscore[:, 2] * zscale), np.min(loadings[2])])
            # ax.set_xlim(min(loadings[0])-0.1, max(loadings[0])+0.1)
            # ax.set_ylim(min(loadings[1])-0.1, max(loadings[1])+0.1)
            # ax.set_zlim(min(loadings[2])-0.1, max(loadings[2])+0.1)
            ax.set_xlim(xlimit_min-0.2, xlimit_max+0.2)
            ax.set_ylim(ylimit_min-0.2, ylimit_max+0.2)
            ax.set_zlim(zlimit_min-0.2, zlimit_max+0.2)
            ax.set_xlabel("PC1 ({}%)".format(var1), fontsize=axlabelfontsize, fontname=axlabelfontname)
            ax.set_ylabel("PC2 ({}%)".format(var2), fontsize=axlabelfontsize, fontname=axlabelfontname)
            ax.set_zlabel("PC3 ({}%)".format(var3), fontsize=axlabelfontsize, fontname=axlabelfontname)
            general.get_figure(show, r, figtype, 'biplot_3d', theme)

    def tsneplot(score=None, axlabelfontsize=9, axlabelfontname="Arial", figtype='png', r=300, show=False,
             markerdot="o", dotsize=6, valphadot=1, colordot='#4a4e4d', colorlist=None, legendpos='best',
             figname='tsne_2d', dim=(6, 4), legendanchor=None, theme=None):
        assert score is not None, "score are missing"
        if theme == 'dark':
            general.dark_bg()
        plt.subplots(figsize=dim)
        if colorlist is not None:
            unique_class = set(colorlist)
            # color_dict = dict()
            assign_values = {col: i for i, col in enumerate(unique_class)}
            color_result_num = [assign_values[i] for i in colorlist]
            if colordot and isinstance(colordot, (tuple, list)):
                colour_map = ListedColormap(colordot)
                s = plt.scatter(score[:, 0], score[:, 1], c=color_result_num, cmap=colour_map,
                                s=dotsize, alpha=valphadot, marker=markerdot)
                plt.legend(handles=s.legend_elements()[0], labels=list(unique_class), loc=legendpos,
                               bbox_to_anchor=legendanchor)
            elif colordot and not isinstance(colordot, (tuple, list)):
                s = plt.scatter(score[:, 0], score[:, 1], c=color_result_num,
                                s=dotsize, alpha=valphadot, marker=markerdot)
                plt.legend(handles=s.legend_elements()[0], labels=list(unique_class), loc=legendpos,
                           bbox_to_anchor=legendanchor)
        else:
            plt.scatter(score[:, 0], score[:, 1], color=colordot,
                       s=dotsize, alpha=valphadot, marker=markerdot)
        plt.xlabel("t-SNE-1", fontsize=axlabelfontsize, fontname=axlabelfontname)
        plt.ylabel("t-SNE-2", fontsize=axlabelfontsize, fontname=axlabelfontname)
        general.get_figure(show, r, figtype, figname, theme)


