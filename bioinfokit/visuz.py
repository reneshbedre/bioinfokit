import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.cm as cmc
import seaborn as sns
from matplotlib_venn import venn3, venn2

def geneplot(d, geneid, lfc, lfc_thr, pv_thr, genenames, gfont, pv):
    if genenames is not None and genenames == "deg":
        for i in d[geneid].unique():
            if (d.loc[d[geneid] == i, lfc].iloc[0] >= lfc_thr and d.loc[d[geneid] == i, pv].iloc[0] < pv_thr) or \
                    (d.loc[d[geneid] == i, lfc].iloc[0] <= -lfc_thr and d.loc[d[geneid] == i, pv].iloc[0] < pv_thr):
                plt.text(d.loc[d[geneid] == i, lfc].iloc[0], d.loc[d[geneid] == i, 'logpv'].iloc[0], i, fontsize=gfont)
    elif genenames is not None and type(genenames) is tuple:
        for i in d[geneid].unique():
            if i in genenames:
                plt.text(d.loc[d[geneid] == i, lfc].iloc[0], d.loc[d[geneid] == i, 'logpv'].iloc[0], i, fontsize=gfont)
    elif genenames is not None and type(genenames) is dict:
        for i in d[geneid].unique():
            if i in genenames:
                plt.text(d.loc[d[geneid] == i, lfc].iloc[0], d.loc[d[geneid] == i, 'logpv'].iloc[0], genenames[i],
                         fontsize=gfont)

def volcano(table="dataset_file", lfc="logFC", pv="p_values", lfc_thr=1, pv_thr=0.05, color=("green", "red"), valpha=1,
            geneid=None, genenames=None, gfont=8):
    # load csv data file
    d = pd.read_csv(table, sep=",")
    color = color
    d.loc[(d[lfc] >= lfc_thr) & (d[pv] < pv_thr), 'color'] = color[0]  # upregulated
    d.loc[(d[lfc] <= -lfc_thr) & (d[pv] < pv_thr), 'color'] = color[1]  # downregulated
    d['color'].fillna('grey', inplace=True)  # intermediate
    d['logpv'] = -(np.log10(d[pv]))

    geneplot(d, geneid, lfc, lfc_thr, pv_thr, genenames, gfont, pv)
    # plot
    plt.scatter(d[lfc], d['logpv'], c=d['color'], alpha=valpha)
    plt.xlabel('log2 Fold Change', fontsize=12, fontname="sans-serif", fontweight="bold")
    plt.ylabel('-log10(P-value)', fontsize=12, fontname="sans-serif", fontweight="bold")
    plt.xticks(fontsize=12, fontname="sans-serif")
    plt.yticks(fontsize=12, fontname="sans-serif")
    
    plt.savefig('volcano.png', format='png', bbox_inches='tight', dpi=300)
    plt.close()


def involcano(table="dataset_file", lfc="logFC", pv="p_values", lfc_thr=1, pv_thr=0.05, color=("green", "red"), valpha=1,
              geneid=None, genenames=None, gfont=8):
    # load csv data file
    d = pd.read_csv(table, sep=",")
    color = color
    d.loc[(d[lfc] >= lfc_thr) & (d[pv] < pv_thr), 'color'] = color[0]  # upregulated
    d.loc[(d[lfc] <= -lfc_thr) & (d[pv] < pv_thr), 'color'] = color[1]  # downregulated
    d['color'].fillna('grey', inplace=True)  # intermediate
    d['logpv'] = -(np.log10(d[pv]))

    geneplot(d, geneid, lfc, lfc_thr, pv_thr, genenames, gfont, pv)

    # plot
    plt.scatter(d[lfc], d['logpv'], c=d['color'], alpha=valpha)
    plt.gca().invert_yaxis()
    plt.xlabel('log2 Fold Change', fontsize=12, fontname="sans-serif", fontweight="bold")
    plt.ylabel('-log10(P-value)', fontsize=12, fontname="sans-serif", fontweight="bold")
    plt.xticks(fontsize=12, fontname="sans-serif")
    plt.yticks(fontsize=12, fontname="sans-serif")
    plt.savefig('involcano.png', format='png', bbox_inches='tight', dpi=300)
    plt.close()


def ma(table="dataset_file", lfc="logFC", ct_count="value1", st_count="value2", lfc_thr=1):
    # load csv data file
    d = pd.read_csv(table, sep=",")
    d.loc[(d[lfc] >= lfc_thr), 'color'] = "green" # upregulated
    d.loc[(d[lfc] <= -lfc_thr), 'color'] = "red"  # downregulated
    d['color'].fillna('grey', inplace=True)  # intermediate
    d['A'] = np.log2((d[ct_count] + d[st_count]) / 2)
    # plot
    plt.scatter(d['A'], d[lfc], c=d['color'])
    # draw a central line at M=0
    plt.axhline(y=0, color='b', linestyle='--')
    plt.xlabel('A', fontsize=15, fontname="sans-serif", fontweight="bold")
    plt.ylabel('M', fontsize=15, fontname="sans-serif", fontweight="bold")
    plt.xticks(fontsize=12, fontname="sans-serif")
    plt.yticks(fontsize=12, fontname="sans-serif")
    plt.savefig('ma.png', format='png', bbox_inches='tight', dpi=300)
    plt.close()


def corr_mat(table="p_df", corm="pearson"):
    # load csv data file
    # d = pd.read_csv(table, sep=",")
    d = pd.DataFrame(data=table)
    d_corr = d.corr(method=corm)
    plt.matshow(d_corr, vmin=-1, vmax=1, cmap=cmc.seismic)
    plt.colorbar()
    cols = list(d)
    ticks = list(range(0, len(list(d))))
    plt.xticks(ticks, cols, fontsize=7, rotation=90)
    plt.yticks(ticks, cols, fontsize=7)
    plt.savefig('corr_mat.png', format='png', bbox_inches='tight', dpi=300)
    plt.close()


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
    # load csv data file
    d = pd.read_csv(table, sep=",")
    d = d.set_index(d.columns[0])
    # plot heatmap without cluster
    # more cmap: https://matplotlib.org/3.1.0/tutorials/colors/colormaps.html
    dim = dim
    fig, hm = plt.subplots(figsize=dim)
    if clus:
        hm = sns.clustermap(d, cmap=cmap, cbar=scale, z_score=zscore, xticklabels=xlabel, yticklabels=ylabel,
                            figsize=dim)
        hm.ax_heatmap.set_xticklabels(hm.ax_heatmap.get_xmajorticklabels(), fontsize=tickfont[0])
        hm.ax_heatmap.set_yticklabels(hm.ax_heatmap.get_ymajorticklabels(), fontsize=tickfont[1])
        plt.savefig('heatmap_clus.png', format='png', bbox_inches='tight', dpi=300)
        plt.close()
    else:
        hm = sns.heatmap(d, cmap=cmap, cbar=scale, xticklabels=xlabel, yticklabels=ylabel)
        plt.xticks(fontsize=tickfont[0])
        plt.yticks(fontsize=tickfont[1])
        plt.savefig('heatmap.png', format='png', bbox_inches='tight', dpi=300)
        plt.close()

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

