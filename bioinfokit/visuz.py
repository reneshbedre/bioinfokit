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
from adjustText import adjust_text


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
                sign_line=False, gstyle=1, show=False):
        # load csv data file
        # d = pd.read_csv(table, sep=",")
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
        general.axis_labels(r'$ -log_{2}(Fold Change)$', r'$ -log_{10}(P-value)$')
        general.get_figure(show, r, 'volcano.png')

    def involcano(d="dataframe", lfc="logFC", pv="p_values", lfc_thr=1, pv_thr=0.05, color=("green", "red"),
                  valpha=1, geneid=None, genenames=None, gfont=8, dim=(5, 5), r=300, ar=90, dotsize=8, markerdot="o",
                sign_line=False, gstyle=1, show=False):
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
        general.axis_labels(r'$ -log_{2}(Fold Change)$', r'$ -log_{10}(P-value)$')
        general.get_figure(show, r, 'involcano.png')

    def ma(df="dataframe", lfc="logFC", ct_count="value1", st_count="value2", lfc_thr=1, valpha=1, dotsize=8,
           markerdot="o", dim=(6, 5), r=300, show=False, color=("green", "red")):
        df.loc[(df[lfc] >= lfc_thr), 'color'] = color[0]  # upregulated
        df.loc[(df[lfc] <= -lfc_thr), 'color'] = color[1]  # downregulated
        df['color'].fillna('grey', inplace=True)  # intermediate
        df['A'] = np.log2((df[ct_count] + df[st_count]) / 2)
        # plot
        plt.subplots(figsize=dim)
        plt.scatter(df['A'], df[lfc], c=df['color'], alpha=valpha, s=dotsize, marker=markerdot)
        # draw a central line at M=0
        plt.axhline(y=0, color='#7d7d7d', linestyle='--')
        general.axis_labels('A', 'M')
        general.get_figure(show, r, 'ma.png')

    def hmap(df="dataframe", cmap="seismic", scale=True, dim=(4, 6), clus=True, zscore=None, xlabel=True,
             ylabel=True, tickfont=(10, 10), r=300, show=False):
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
            general.get_figure(show, r, 'heatmap.png')
        else:
            hm = sns.heatmap(df, cmap=cmap, cbar=scale, xticklabels=xlabel, yticklabels=ylabel)
            plt.xticks(fontsize=tickfont[0])
            plt.yticks(fontsize=tickfont[1])
            general.get_figure(show, r, 'heatmap.png')


class general():
    def __init__(self):
        pass

    rand_colors = ('#a7414a', '#282726', '#6a8a82', '#a37c27', '#563838', '#0584f2', '#f28a30', '#f05837',
                   '#6465a5', '#00743f', '#be9063', '#de8cf0', '#888c46', '#c0334d', '#270101', '#8d2f23',
                   '#ee6c81', '#65734b', '#14325c', '#704307', '#b5b3be', '#f67280', '#ffd082', '#ffd800',
                   '#ad62aa', '#21bf73', '#a0855b', '#5edfff', '#08ffc8', '#ca3e47', '#c9753d', '#6c5ce7')

    def get_figure(show, r, fig_name):
        if show:
            plt.show()
        else:
            plt.savefig(fig_name, format='png', bbox_inches='tight', dpi=r)
        plt.close()

    def axis_labels(x, y):
        plt.xlabel(x, fontsize=12, fontname="sans-serif")
        plt.ylabel(y, fontsize=12, fontname="sans-serif")
        plt.xticks(fontsize=9, fontname="sans-serif")
        plt.yticks(fontsize=9, fontname="sans-serif")

    def depr_mes(func_name):
        print("This function is deprecated. Please use", func_name )
        print("Read docs at https://reneshbedre.github.io/blog/howtoinstall.html")

class marker():
    def geneplot_mhat(df, markeridcol, chr, pv, gwasp, markernames, gfont, ax):
        if markeridcol is not None:
            if markernames is not None and markernames is True:
                for i in df[markeridcol].unique():
                    if df.loc[df[markeridcol] == i, pv].iloc[0] <= gwasp:
                        texts = [plt.text((df.loc[df[markeridcol] == i, 'ind'].iloc[0]), df.loc[df[markeridcol] == i, 'tpval'].iloc[0],
                                 str(i), fontsize=gfont)]
                        adjust_text(texts)
            elif markernames is not None and type(markernames) is tuple:
                for i in df[markeridcol].unique():
                    if i in markernames:
                        texts = [plt.text(df.loc[df[markeridcol] == i, 'ind'].iloc[0], df.loc[df[markeridcol] == i, 'tpval'].iloc[0],
                                 str(i), fontsize=gfont)]
                        adjust_text(texts)
            elif markernames is not None and type(markernames) is dict:
                for i in df[markeridcol].unique():
                    if i in markernames:
                        texts = [plt.text(df.loc[df[markeridcol] == i, 'ind'].iloc[0], df.loc[df[markeridcol] == i, 'tpval'].iloc[0],
                                 markernames[i], fontsize=gfont)]
                        adjust_text(texts)
        else:
            print("Error: provide 'markeridcol' parameter")
            sys.exit(1)

    def mhat(df="dataframe", chr=None, pv=None, color=None, dim=(6,4), r=300, ar=90, gwas_sign_line=False,
             gwasp=5E-08, dotsize=8, markeridcol=None, markernames=None, gfont=8, valpha=1, show=False):

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
        ax.set_yticks(np.arange(0, max(df['tpval']+1), 1))
        ax.set_xticklabels(xlabels, fontsize=9, rotation=ar)
        ax.set_xlabel('Chromosomes', fontsize=9, fontname="sans-serif")
        ax.set_ylabel(r'$\bf -log_{10}(P)$', fontsize=9, fontname="sans-serif", fontweight="bold")
        ax.set_ylim([0, max(df['tpval']+1)])
        general.get_figure(show, r, 'manhatten.png')


class stat:
    def __init__(self):
        pass

    def bardot(df="dataframe", dim=(6, 4), bw=0.4, colorbar="#f2aa4cff", colordot=["#101820ff"], hbsize=4, r=300, ar=0,
               dotsize=6, valphabar=1, valphadot=1, markerdot="o", errorbar=True, show=False):
        xbar = np.arange(len(df.columns.to_numpy()))
        color_list_bar = colorbar
        color_list_dot = colordot
        if len(color_list_dot) == 1:
            color_list_dot = colordot*len(df.columns.to_numpy())
        fig, ax = plt.subplots(figsize=dim)
        if errorbar:
            ax.bar(x=xbar, height=df.describe().loc['mean'], yerr=df.sem(), width=bw, color=color_list_bar, capsize=hbsize,
                zorder=0, alpha=valphabar)
        else:
            ax.bar(x=xbar, height=df.describe().loc['mean'], width=bw, color=color_list_bar,
                   capsize=hbsize,
                   zorder=0, alpha=valphabar)
        ax.set_xticks(xbar)
        ax.set_xticklabels(df.columns.to_numpy(), fontsize=9, rotation=ar)
        # add dots
        for cols in range(len(df.columns.to_numpy())):
            # get markers from here https://matplotlib.org/3.1.1/api/markers_api.html
            ax.scatter(x=np.linspace(xbar[cols]-bw/2, xbar[cols]+bw/2, int(df.describe().loc['count'][cols])),
                       y=df[df.columns[cols]], s=dotsize, color=color_list_dot[cols], zorder=1, alpha=valphadot,
                       marker=markerdot)
        general.get_figure(show, r, 'bardot.png')

    def regplot(df="dataframe", x=None, y=None, yhat=None, dim=(6, 4), colordot='#4a4e4d', colorline='#fe8a71', r=300,
                ar=0, dotsize=6, valphaline=1, valphadot=1, linewidth=1, markerdot="o", show=False):
        fig, ax = plt.subplots(figsize=dim)
        plt.scatter(df[x].to_numpy(), df[y].to_numpy(), color=colordot, s=dotsize, alpha=valphadot, marker=markerdot,
                    label='Observed data')
        plt.plot(df[x].to_numpy(), df[yhat].to_numpy(), color=colorline, linewidth=linewidth, alpha=valphaline,
                 label='Regression line')
        plt.xlabel(x, fontsize=9)
        plt.ylabel(y, fontsize=9)
        plt.legend(fontsize=8)
        general.get_figure(show, r, 'reg_plot.png')

    def reg_resid_plot(df="dataframe", yhat=None, resid=None, stdresid=None, dim=(6, 4), colordot='#4a4e4d',
                       colorline='#2ab7ca', r=300, ar=0, dotsize=6, valphaline=1, valphadot=1, linewidth=1,
                       markerdot="o", show=False):
        fig, ax = plt.subplots(figsize=dim)
        if resid is not None:
            plt.scatter(df[yhat], df[resid], color=colordot, s=dotsize, alpha=valphadot, marker=markerdot)
            plt.axhline(y=0, color=colorline, linestyle='--', linewidth=linewidth, alpha=valphaline)
            plt.xlabel("Fitted")
            plt.ylabel("Residuals")
            general.get_figure(show, r, 'resid_plot.png')
        else:
            print ("Error: Provide residual data")
        if stdresid is not None:
            plt.scatter(df[yhat], df[stdresid], color=colordot, s=dotsize, alpha=valphadot, marker=markerdot)
            plt.axhline(y=0, color=colorline, linestyle='--', linewidth=linewidth, alpha=valphaline)
            plt.xlabel("Fitted")
            plt.ylabel("Standardized Residuals")
            general.get_figure(show, r, 'std_resid_plot.png')
        else:
            print ("Error: Provide standardized residual data")

    def corr_mat(df="dataframe", corm="pearson", cmap="seismic", r=300, show=False, dim=(6, 5)):
        d_corr = df.corr(method=corm)
        plt.subplots(figsize=dim)
        plt.matshow(d_corr, vmin=-1, vmax=1, cmap=cmap)
        plt.colorbar()
        cols = list(df)
        ticks = list(range(0, len(list(df))))
        plt.xticks(ticks, cols, fontsize=7, rotation=90)
        plt.yticks(ticks, cols, fontsize=7)
        general.get_figure(show, r, 'corr_mat.png')

class help:
    def __init__(self):
        pass

    @staticmethod
    def mhat():
        text = """
        Manhatten plot

        bioinfokit.visuz.marker.mhat(df, chr, pv, color, dim, r, ar, gwas_sign_line, gwasp, dotsize, markeridcol, markernames, gfont, valpha)

        Parameters:
        ------------
        df             : Pandas dataframe object with atleast SNP, chromosome, and P-values columns
        chr            : Name of a column having chromosome numbers [string][default:None]
        pv             : Name of a column having P-values. Must be numeric column [string][default:None]
        color          : List the name of the colors to be plotted. It can accept two alternate colors or the number colors 
                         equal to chromosome number. If nothing (None) provided, it will randomly assign the color to each 
                         chromosome [list][default:None]
        dim            : Figure size [tuple of two floats (width, height) in inches][default: (6, 4)]
        r              : Figure resolution in dpi [int][default: 300]
        ar             : Rotation of X-axis labels [float][default: 90]
        gwas_sign_line : Plot statistical significant threshold line defined by option `gwasp` 
                         [bool (True or False)][default: False]
        gwasp          : Statistical significant threshold to identify significant SNPs [float][default: 5E-08]
        dotsize        : The size of the dots in the plot [float][default: 8]
        markeridcol    : Name of a column having SNPs. This is necessary for plotting SNP names on the plot 
                         [string][default: None]
        markernames    : The list of the SNPs to display on the plot. These SNP should be present in SNP column. 
                         Additionally, it also accepts the dict of SNPs and its associated gene name. If this option set 
                         to True, it will label all SNPs with P-value significant score defined by `gwasp` 
                         [string, list, dict][default: True]
        gfont          : Font size for SNP names to display on the plot [float][default: 8]
        valpha         : Transparency of points on plot [float (between 0 and 1)][default: 1.0]

        Returns:
        Manhatten plot image in same directory (manhatten.png)

        Working example: https://reneshbedre.github.io/blog/manhat.html
        """

        print(text)

    def bardot():
        text = """
               Bar-dot plot

               bioinfokit.visuz.stat.bardot(df, colorbar, colordot, bw, dim, r, ar, hbsize, errorbar, dotsize, markerdot, valphabar, valphadot)

               Parameters:
               ------------
               df        : Pandas dataframe object
               colorbar  : Color of bar graph [string or list][default:"#bbcfff"]
               colordot  : Color of dots on bar [string or list][default:"#ee8972"]
               bw        : Width of bar [float][default: 0.4]
               dim       : Figure size [tuple of two floats (width, height) in inches][default: (6, 4)]
               r         : Figure resolution in dpi [int][default: 300]
               ar        : Rotation of X-axis labels [float][default: 0]
               hbsize    : Horizontal bar size for standard error bars [float][default: 4]
               errorbar  : Draw standard error bars [bool (True or False)][default: True]
               dotsize   : The size of the dots in the plot [float][default: 6]
               markerdot : Shape of the dot marker. 
                           See more options at  https://matplotlib.org/3.1.1/api/markers_api.html [string][default: "o"]
               valphabar : Transparency of bars on plot [float (between 0 and 1)][default: 1]
               valphadot : Transparency of dots on plot [float (between 0 and 1)][default: 1]

               Returns:
               Bar-dot plot image in same directory (bardot.png)

               Working example: https://reneshbedre.github.io/blog/bardot.html
               """

        print(text)