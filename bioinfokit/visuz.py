import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cmc


def volcano(table="dataset_file", lfc="logFC", pv="p_values", lfc_thr=1, pv_thr=0.05):
    # load csv data file
    d = pd.read_csv(table, sep=",")
    d.loc[(d[lfc] >= lfc_thr) & (d[pv] < pv_thr), 'color'] = "green"  # upregulated
    d.loc[(d[lfc] <= -lfc_thr) & (d[pv] < pv_thr), 'color'] = "red"  # downregulated
    d['color'].fillna('grey', inplace=True)  # intermediate
    d['logpv'] = -(np.log10(d[pv]))
    # plot
    plt.scatter(d[lfc], d['logpv'], c=d['color'])
    plt.xlabel('log2 Fold Change', fontsize=15, fontname="sans-serif", fontweight="bold")
    plt.ylabel('-log10(P-value)', fontsize=15, fontname="sans-serif", fontweight="bold")
    plt.xticks(fontsize=12, fontname="sans-serif")
    plt.yticks(fontsize=12, fontname="sans-serif")
    plt.savefig('volcano.png', format='png', bbox_inches='tight', dpi=300)
    plt.close()


def involcano(table="dataset_file", lfc="logFC", pv="p_values", lfc_thr=1, pv_thr=0.05):
    # load csv data file
    d = pd.read_csv(table, sep=",")
    d.loc[(d[lfc] >= lfc_thr) & (d[pv] < pv_thr), 'color'] = "green"  # upregulated
    d.loc[(d[lfc] <= -lfc_thr) & (d[pv] < pv_thr), 'color'] = "red"  # downregulated
    d['color'].fillna('grey', inplace=True)  # intermediate
    d['logpv'] = -(np.log10(d[pv]))
    # plot
    plt.scatter(d[lfc], d['logpv'], c=d['color'])
    plt.gca().invert_yaxis()
    plt.xlabel('log2 Fold Change', fontsize=15, fontname="sans-serif", fontweight="bold")
    plt.ylabel('-log10(P-value)', fontsize=15, fontname="sans-serif", fontweight="bold")
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


