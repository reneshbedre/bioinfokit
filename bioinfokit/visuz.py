import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.cm as cmc
import seaborn as sns
from matplotlib_venn import venn3, venn2
from random import sample
from functools import reduce
import sys



def volcano(d="dataframe", lfc=None, pv=None, lfc_thr=1, pv_thr=0.05, color=("green", "red"), valpha=1,
            geneid=None, genenames=None, gfont=8, dim=(6,4), r=300, ar=90, dotsize=8, markerdot="o", sign_line=False):
    general.depr_mes("bioinfokit.visuz.gene_exp.volcano")


def involcano(table="dataset_file", lfc="logFC", pv="p_values", lfc_thr=1, pv_thr=0.05, color=("green", "red"), valpha=1,
              geneid=None, genenames=None, gfont=8):
    general.depr_mes("bioinfokit.visuz.gene_exp.involcano")


def ma(table="dataset_file", lfc="logFC", ct_count="value1", st_count="value2", lfc_thr=1):
    general.depr_mes("bioinfokit.visuz.gene_exp.ma")

def corr_mat(table="p_df", corm="pearson"):
    general.depr_mes("bioinfokit.visuz.stat.corr_mat")

def screeplot(obj="pcascree"):
    y = [x * 100 for x in obj[1]]
    plt.bar(obj[0], y)
    plt.xlabel('PCs', fontsize=12, fontname="sans-serif")
    plt.ylabel('Proportion of variance (%)', fontsize=12, fontname="sans-serif")
    plt.xticks(fontsize=7, rotation=70)
    plt.savefig('screeplot.png', format='png', bbox_inches='tight', dpi=300)
    plt.close()


def pcaplot(x="x", y="y", z="z", labels="d_cols", var1="var1", var2="var2", var3="var3"):
    for i, varnames in enumerate(labels):
        plt.scatter(x[i], y[i])
        plt.text(x[i], y[i], varnames, fontsize=10)
    plt.xlabel("PC1 ({}%)".format(var1), fontsize=12, fontname="sans-serif")
    plt.ylabel("PC2 ({}%)".format(var2), fontsize=12, fontname="sans-serif")
    plt.tight_layout()
    plt.savefig('pcaplot_2d.png', format='png', bbox_inches='tight', dpi=300)
    plt.close()

    # for 3d plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    for i, varnames in enumerate(labels):
        ax.scatter(x[i], y[i], z[i])
        ax.text(x[i], y[i], z[i], varnames, fontsize=10)
    ax.set_xlabel("PC1 ({}%)".format(var1), fontsize=12, fontname="sans-serif")
    ax.set_ylabel("PC2 ({}%)".format(var2), fontsize=12, fontname="sans-serif")
    ax.set_zlabel("PC3 ({}%)".format(var3), fontsize=12, fontname="sans-serif")
    plt.tight_layout()
    plt.savefig('pcaplot_3d.png', format='png', bbox_inches='tight',  dpi=300)
    plt.close()


def hmap(table="dataset_file", cmap="seismic", scale=True, dim=(4,6), clus=True, zscore=None, xlabel=True, ylabel=True,
         tickfont=(10,10)):
    general.depr_mes("bioinfokit.visuz.gene_exp.hmap")

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

class gene_exp():
    def __init__(self):
        pass

    def geneplot(d, geneid, lfc, lfc_thr, pv_thr, genenames, gfont, pv, gstyle):
        if genenames is not None and genenames == "deg":
            for i in d[geneid].unique():
                if (d.loc[d[geneid] == i, lfc].iloc[0] >= lfc_thr and d.loc[d[geneid] == i, pv].iloc[0] < pv_thr) or \
                        (d.loc[d[geneid] == i, lfc].iloc[0] <= -lfc_thr and d.loc[d[geneid] == i, pv].iloc[0] < pv_thr):
                    if gstyle==1:
                        plt.text(d.loc[d[geneid] == i, lfc].iloc[0], d.loc[d[geneid] == i, 'logpv'].iloc[0], i,
                                      fontsize=gfont)
                    elif gstyle==2:
                        plt.annotate(i, xy=(d.loc[d[geneid] == i, lfc].iloc[0], d.loc[d[geneid] == i, 'logpv'].iloc[0]),
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
                        plt.text(d.loc[d[geneid] == i, lfc].iloc[0], d.loc[d[geneid] == i, 'logpv'].iloc[0], i,
                                      fontsize=gfont)
                    elif gstyle==2:
                        plt.annotate(i, xy=(d.loc[d[geneid] == i, lfc].iloc[0], d.loc[d[geneid] == i, 'logpv'].iloc[0]),
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
                        plt.text(d.loc[d[geneid] == i, lfc].iloc[0], d.loc[d[geneid] == i, 'logpv'].iloc[0],
                                      genenames[i], fontsize=gfont)
                    elif gstyle == 2:
                        plt.annotate(genenames[i], xy=(d.loc[d[geneid] == i, lfc].iloc[0], d.loc[d[geneid] == i, 'logpv'].iloc[0]),
                                     xycoords='data', xytext=(5, -15), textcoords='offset points', size=6,
                                     bbox=dict(boxstyle="round", alpha=0.1),
                                     arrowprops=dict(arrowstyle="wedge,tail_width=0.5", alpha=0.1, relpos=(0, 0)))
                    else:
                        print("Error: invalid gstyle choice")
                        sys.exit(1)

    def volcano(d="dataframe", lfc=None, pv=None, lfc_thr=1, pv_thr=0.05, color=("green", "red"), valpha=1,
                geneid=None, genenames=None, gfont=8, dim=(5, 5), r=300, ar=90, dotsize=8, markerdot="o",
                sign_line=False, gstyle=1, show=False, figtype='png', axtickfontsize=9,
               axtickfontname="Arial", axlabelfontsize=9, axlabelfontname="Arial", axxlabel=None,
                axylabel=None, xlm=None, ylm=None):
        _x = r'$ log_{2}(Fold Change)$'
        _y = r'$ -log_{10}(P-value)$'
        color = color
        d.loc[(d[lfc] >= lfc_thr) & (d[pv] < pv_thr), 'color'] = color[0]  # upregulated
        d.loc[(d[lfc] <= -lfc_thr) & (d[pv] < pv_thr), 'color'] = color[1]  # downregulated
        d['color'].fillna('grey', inplace=True)  # intermediate
        d['logpv'] = -(np.log10(d[pv]))
        # plot
        plt.subplots(figsize=dim)
        plt.scatter(d[lfc], d['logpv'], c=d['color'], alpha=valpha, s=dotsize, marker=markerdot)
        if sign_line:
            plt.axhline(y=-np.log10(pv_thr), linestyle='--', color='#7d7d7d', linewidth=1)
            plt.axvline(x=lfc_thr, linestyle='--', color='#7d7d7d', linewidth=1)
            plt.axvline(x=-lfc_thr, linestyle='--', color='#7d7d7d', linewidth=1)
        gene_exp.geneplot(d, geneid, lfc, lfc_thr, pv_thr, genenames, gfont, pv, gstyle)
        if axxlabel:
            _x = axxlabel
        if axylabel:
            _y = axylabel
        general.axis_labels(_x, _y, axlabelfontsize, axlabelfontname)
        general.axis_ticks(xlm, ylm, axtickfontsize, axtickfontname, ar)
        general.get_figure(show, r, figtype, 'volcano')

    def involcano(d="dataframe", lfc="logFC", pv="p_values", lfc_thr=1, pv_thr=0.05, color=("green", "red"),
                  valpha=1, geneid=None, genenames=None, gfont=8, dim=(5, 5), r=300, ar=90, dotsize=8, markerdot="o",
                sign_line=False, gstyle=1, show=False, figtype='png', axtickfontsize=9,
               axtickfontname="Arial", axlabelfontsize=9, axlabelfontname="Arial", axxlabel=None,
                axylabel=None, xlm=None, ylm=None):
        _x = r'$ log_{2}(Fold Change)$'
        _y = r'$ -log_{10}(P-value)$'
        color = color
        d.loc[(d[lfc] >= lfc_thr) & (d[pv] < pv_thr), 'color'] = color[0]  # upregulated
        d.loc[(d[lfc] <= -lfc_thr) & (d[pv] < pv_thr), 'color'] = color[1]  # downregulated
        d['color'].fillna('grey', inplace=True)  # intermediate
        d['logpv'] = -(np.log10(d[pv]))

        # plot
        plt.subplots(figsize=dim)
        plt.scatter(d[lfc], d['logpv'], c=d['color'], alpha=valpha, s=dotsize, marker=markerdot)
        gene_exp.geneplot(d, geneid, lfc, lfc_thr, pv_thr, genenames, gfont, pv, gstyle)
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
        general.get_figure(show, r, figtype, 'involcano')

    def ma(df="dataframe", lfc="logFC", ct_count="value1", st_count="value2", lfc_thr=1, valpha=1, dotsize=8,
           markerdot="o", dim=(6, 5), r=300, show=False, color=("green", "red"), ar=90, figtype='png', axtickfontsize=9,
           axtickfontname="Arial", axlabelfontsize=9, axlabelfontname="Arial", axxlabel=None,
           axylabel=None, xlm=None, ylm=None):
        _x, _y = 'A', 'M'
        df.loc[(df[lfc] >= lfc_thr), 'color'] = color[0]  # upregulated
        df.loc[(df[lfc] <= -lfc_thr), 'color'] = color[1]  # downregulated
        df['color'].fillna('grey', inplace=True)  # intermediate
        df['A'] = np.log2((df[ct_count] + df[st_count]) / 2)
        # plot
        plt.subplots(figsize=dim)
        plt.scatter(df['A'], df[lfc], c=df['color'], alpha=valpha, s=dotsize, marker=markerdot)
        # draw a central line at M=0
        plt.axhline(y=0, color='#7d7d7d', linestyle='--')
        if axxlabel:
            _x = axxlabel
        if axylabel:
            _y = axylabel
        general.axis_labels(_x, _y, axlabelfontsize, axlabelfontname)
        general.axis_ticks(xlm, ylm, axtickfontsize, axtickfontname, ar)
        general.get_figure(show, r, figtype, 'ma')

    def hmap(df="dataframe", cmap="seismic", scale=True, dim=(4, 6), clus=True, zscore=None, xlabel=True,
             ylabel=True, tickfont=(10, 10), r=300, show=False, figtype='png'):
        # df = df.set_index(d.columns[0])
        # plot heatmap without cluster
        # more cmap: https://matplotlib.org/3.1.0/tutorials/colors/colormaps.html
        # dim = dim
        fig, hm = plt.subplots(figsize=dim)
        if clus:
            hm = sns.clustermap(df, cmap=cmap, cbar=scale, z_score=zscore, xticklabels=xlabel, yticklabels=ylabel,
                                figsize=dim)
            hm.ax_heatmap.set_xticklabels(hm.ax_heatmap.get_xmajorticklabels(), fontsize=tickfont[0])
            hm.ax_heatmap.set_yticklabels(hm.ax_heatmap.get_ymajorticklabels(), fontsize=tickfont[1])
            general.get_figure(show, r, figtype, 'heatmap')
        else:
            hm = sns.heatmap(df, cmap=cmap, cbar=scale, xticklabels=xlabel, yticklabels=ylabel)
            plt.xticks(fontsize=tickfont[0])
            plt.yticks(fontsize=tickfont[1])
            general.get_figure(show, r, figtype, 'heatmap')


class general():
    def __init__(self):
        pass

    rand_colors = ('#a7414a', '#282726', '#6a8a82', '#a37c27', '#563838', '#0584f2', '#f28a30', '#f05837',
                   '#6465a5', '#00743f', '#be9063', '#de8cf0', '#888c46', '#c0334d', '#270101', '#8d2f23',
                   '#ee6c81', '#65734b', '#14325c', '#704307', '#b5b3be', '#f67280', '#ffd082', '#ffd800',
                   '#ad62aa', '#21bf73', '#a0855b', '#5edfff', '#08ffc8', '#ca3e47', '#c9753d', '#6c5ce7')

    def get_figure(show, r, figtype, fig_name):
        if show:
            plt.show()
        else:
            plt.savefig(fig_name+'.'+figtype, format=figtype, bbox_inches='tight', dpi=r)
        plt.close()

    def axis_labels(x, y, axlabelfontsize=None, axlabelfontname=None):
        plt.xlabel(x, fontsize=axlabelfontsize, fontname=axlabelfontname)
        plt.ylabel(y, fontsize=axlabelfontsize, fontname=axlabelfontname)
        # plt.xticks(fontsize=9, fontname="sans-serif")
        # plt.yticks(fontsize=9, fontname="sans-serif")

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

    def depr_mes(func_name):
        print("This function is deprecated. Please use", func_name )
        print("Read docs at https://reneshbedre.github.io/blog/howtoinstall.html")

class marker():
    def geneplot_mhat(df, markeridcol, chr, pv, gwasp, markernames, gfont, ax):
        if markeridcol is not None:
            if markernames is not None and markernames is True:
                for i in df[markeridcol].unique():
                    if df.loc[df[markeridcol] == i, pv].iloc[0] <= gwasp:
                        plt.text((df.loc[df[markeridcol] == i, 'ind'].iloc[0]), df.loc[df[markeridcol] == i, 'tpval'].iloc[0],
                                 str(i), fontsize=gfont)
            elif markernames is not None and type(markernames) is tuple:
                for i in df[markeridcol].unique():
                    plt.text(df.loc[df[markeridcol] == i, 'ind'].iloc[0], df.loc[df[markeridcol] == i, 'tpval'].iloc[0],
                                 str(i), fontsize=gfont)
            elif markernames is not None and type(markernames) is dict:
                for i in df[markeridcol].unique():
                    if i in markernames:
                       plt.text(df.loc[df[markeridcol] == i, 'ind'].iloc[0], df.loc[df[markeridcol] == i, 'tpval'].iloc[0],
                                 markernames[i], fontsize=gfont)
        else:
            print("Error: provide 'markeridcol' parameter")
            sys.exit(1)

    def mhat(df="dataframe", chr=None, pv=None, color=None, dim=(6,4), r=300, ar=90, gwas_sign_line=False,
             gwasp=5E-08, dotsize=8, markeridcol=None, markernames=None, gfont=8, valpha=1, show=False, figtype='png',
             axxlabel=None, axylabel=None, axlabelfontsize=9, axlabelfontname="Arial", axtickfontsize=9,
             axtickfontname="Arial", ylm=None):

        _x, _y = 'Chromosomes', r'$ -log_{10}(P)$'
        rand_colors = ('#a7414a', '#282726', '#6a8a82', '#a37c27', '#563838', '#0584f2', '#f28a30', '#f05837',
                       '#6465a5', '#00743f', '#be9063', '#de8cf0', '#888c46', '#c0334d', '#270101', '#8d2f23',
                       '#ee6c81', '#65734b', '#14325c', '#704307', '#b5b3be', '#f67280', '#ffd082', '#ffd800',
                       '#ad62aa', '#21bf73', '#a0855b', '#5edfff', '#08ffc8', '#ca3e47', '#c9753d', '#6c5ce7')
        '''
         rand_colors = ('#f67280', '#00a8cc', '#ffd082', '#fb8d62', '#6e5773', '#21bf73', '#d5c455', '#c9753d',
                       '#ad62aa','#d77fa1', '#a0855b', '#ffd800', '#da2d2d', '#6f9a8d', '#a8ff3e', '#b2fcff',
                       '#a0c334', '#b5525c', '#c06c84', '#3a3535', '#9b45e4', '#f6da63', '#9dab86', '#0c093c',
                       '#f6f078', '#64c4ed', '#da4302', '#5edfff', '#08ffc8', '#ca3e47', '#f7ff56', '#6c5ce7')
        '''
        # minus log10 of P-value
        df['tpval'] = -np.log10(df[pv])
        df = df.sort_values(chr)
        # add indices
        df['ind'] = range(len(df))
        df_group = df.groupby(chr)
        if color is not None and len(color) == 2:
            color_1 = int(df[chr].nunique() / 2) * [color[0]]
            color_2 = int(df[chr].nunique() / 2) * [color[1]]
            if df[chr].nunique() % 2 == 0:
                color_list = list(reduce(lambda x, y: x+y, zip(color_1, color_2)))
            elif df[chr].nunique() % 2 == 1:
                color_list = list(reduce(lambda x, y: x + y, zip(color_1, color_2)))
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
            marker.geneplot_mhat(df, markeridcol, chr, pv, gwasp, markernames, gfont, ax=ax)
        ax.margins(x=0)
        ax.margins(y=0)
        ax.set_xticks(xticks)
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
        general.get_figure(show, r, figtype, 'manhatten')


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
                       y=df[df.columns[cols]], s=dotsize, color=color_list_dot[cols], zorder=1, alpha=valphadot,
                       marker=markerdot)
        general.get_figure(show, r, figtype, 'bardot')

    def regplot(df="dataframe", x=None, y=None, yhat=None, dim=(6, 4), colordot='#4a4e4d', colorline='#fe8a71', r=300,
                ar=0, dotsize=6, valphaline=1, valphadot=1, linewidth=1, markerdot="o", show=False, axtickfontsize=9,
               axtickfontname="Arial", axlabelfontsize=9, axlabelfontname="Arial", ylm=None, xlm=None, axxlabel=None,
                axylabel=None, figtype='png'):
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
        general.get_figure(show, r, figtype, 'reg_plot')

    def reg_resid_plot(df="dataframe", yhat=None, resid=None, stdresid=None, dim=(6, 4), colordot='#4a4e4d',
                       colorline='#2ab7ca', r=300, ar=0, dotsize=6, valphaline=1, valphadot=1, linewidth=1,
                       markerdot="o", show=False, figtype='png'):
        fig, ax = plt.subplots(figsize=dim)
        if resid is not None:
            plt.scatter(df[yhat], df[resid], color=colordot, s=dotsize, alpha=valphadot, marker=markerdot)
            plt.axhline(y=0, color=colorline, linestyle='--', linewidth=linewidth, alpha=valphaline)
            plt.xlabel("Fitted")
            plt.ylabel("Residuals")
            general.get_figure(show, r, figtype, 'resid_plot')
        else:
            print ("Error: Provide residual data")
        if stdresid is not None:
            plt.scatter(df[yhat], df[stdresid], color=colordot, s=dotsize, alpha=valphadot, marker=markerdot)
            plt.axhline(y=0, color=colorline, linestyle='--', linewidth=linewidth, alpha=valphaline)
            plt.xlabel("Fitted")
            plt.ylabel("Standardized Residuals")
            general.get_figure(show, r, figtype, 'std_resid_plot')
        else:
            print ("Error: Provide standardized residual data")

    def corr_mat(df="dataframe", corm="pearson", cmap="seismic", r=300, show=False, dim=(6, 5), axtickfontname="Arial",
                 axtickfontsize=7, ar=90, figtype='png'):
        d_corr = df.corr(method=corm)
        plt.subplots(figsize=dim)
        plt.matshow(d_corr, vmin=-1, vmax=1, cmap=cmap)
        plt.colorbar()
        cols = list(df)
        ticks = list(range(0, len(list(df))))
        plt.xticks(ticks, cols, fontsize=axtickfontsize, fontname=axtickfontname, rotation=ar)
        plt.yticks(ticks, cols, fontsize=axtickfontsize, fontname=axtickfontname)
        general.get_figure(show, r, figtype, 'corr_mat')


