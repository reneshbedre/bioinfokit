from sklearn.decomposition import PCA
import pandas as pd
import re
import numpy as np
from bioinfokit.visuz import screeplot, pcaplot
from itertools import groupby
import string
import sys, csv
import matplotlib.pyplot as plt
import scipy.stats as stats
from tabulate import tabulate
from termcolor import colored
from statsmodels.graphics.mosaicplot import mosaic


def seqcov(file="fastq_file", gs="genome_size"):
    x = fastq_format_check(file)
    if x == 1:
        print("Error: Sequences are not in fastq format")
        sys.exit(1)
    num_reads, total_len = fqreadcounter(file)
    # haploid genome_size must be in Mbp; convert in bp
    gs = gs * 1e6
    cov = round(float(total_len / gs), 2)
    print("Sequence coverage for", file, "is", cov)

def mergevcf(file="vcf_file_com_sep"):
    vcf_files = file.split(",")
    merge_vcf = open("merge_vcf.vcf", "w+")
    file_count = 0
    print("merging vcf files...")
    for f in vcf_files:
        if file_count == 0:
            read_file = open(f, "rU")
            for line in read_file:
                merge_vcf.write(line)
            read_file.close()
        elif file_count > 0:
            read_file = open(f, "rU")
            for line in read_file:
                if not line.startswith("#"):
                    merge_vcf.write(line)
            read_file.close()
        file_count += 1
    merge_vcf.close()


def pca(table="p_df"):
    d = pd.DataFrame(data=table)
    d_cols = list(d.columns.values)
    pca_out = PCA()
    pca_out.fit(d)
    prop_var = pca_out.explained_variance_ratio_
    cum_prop_var = np.cumsum(prop_var)
    rotation = pca_out.components_
    num_pc = pca_out.n_features_
    pc_list = list(range(1, num_pc+1))
    pc_list = ["PC"+str(w) for w in pc_list]
    pca_df_var = [prop_var, cum_prop_var]
    pca_df_out = pd.DataFrame.from_dict(dict(zip(pc_list, zip(*pca_df_var))))
    pca_df_rot_out = pd.DataFrame.from_dict(dict(zip(pc_list, rotation)))
    pca_df_out.rename(index={0: "Proportion of Variance", 1: "Cumulative proportion"}, inplace=True)
    print("Component summary\n")
    print(pca_df_out)
    print("\nLoadings\n")
    pca_df_rot_out['sample'] = d_cols
    pca_df_rot_out = pca_df_rot_out.set_index('sample')
    del pca_df_rot_out.index.name
    print(pca_df_rot_out)
    pcascree = [pc_list, prop_var]
    # screeplot
    screeplot(obj=pcascree)
    # for pcaplot; take PC1 and PC2 loadings
    pcaplot(x=rotation[0], y=rotation[1], z=rotation[2], labels=d_cols, var1=round(prop_var[0]*100, 2), var2=round(prop_var[1]*100, 2),
            var3=round(prop_var[2] * 100, 2))


def extract_seq(file="fasta_file", id="id_file"):
    # extract seq from fasta file based on id match
    id_list = []
    id_file = open(id, "rU")
    out_file = open("output.fasta", 'w')
    for line in id_file:
        id_name = line.rstrip('\n')
        id_list.append(id_name)
    list_len = len(id_list)
    value = [1] * list_len
    # id_list converted to dict for faster search
    dict_list = dict(zip(id_list, value))
    fasta_iter = fasta_reader(file)
    for record in fasta_iter:
        fasta_header, seq = record
        if fasta_header.strip() in dict_list.keys():
            out_file.write(">"+fasta_header+"\n"+seq+"\n")
    out_file.close()
    id_file.close()


# remove seqs which match to ids in id file
def extract_seq_nomatch(file="fasta_file", id="id_file"):
    # extract seq from fasta file based on id match
    id_list = []
    id_file = open(id, "rU")
    out_file = open("output.fasta", 'w')
    for line in id_file:
        id_name = line.rstrip('\n')
        id_list.append(id_name)
    list_len = len(id_list)
    value = [1] * list_len
    # id_list converted to dict for faster search
    dict_list = dict(zip(id_list, value))
    fasta_iter = fasta_reader(file)
    for record in fasta_iter:
        fasta_header, seq = record
        if fasta_header.strip() not in dict_list.keys():
            out_file.write(">"+fasta_header+"\n"+seq+"\n")
    out_file.close()
    id_file.close()

def fqreadcounter(file="fastq_file"):
    read_file = open(file, "rU")
    num_lines = 0
    total_len = 0
    for line in read_file:
        num_lines += 1
        header_1 = line.rstrip()
        read = next(read_file).rstrip()
        len_read = len(read)
        total_len += len_read
        header_2 = next(read_file).rstrip()
        read_qual = next(read_file).rstrip()
    read_file.close()
    num_reads = num_lines/4
    return num_reads, total_len


def fasta_reader(file="fasta_file"):
    read_file = open(file, "rU")
    fasta_iter = (rec[1] for rec in groupby(read_file, lambda line: line[0] == ">"))
    for record in fasta_iter:
        fasta_header = record .__next__()[1:].strip()
        fasta_header = re.split("\s+", fasta_header)[0]
        seq = "".join(s.strip() for s in fasta_iter.__next__())
        yield (fasta_header, seq)

def rev_com(seq=None, file=None):
    if seq is not None:
        rev_seq = seq[::-1]
        rev_seq = rev_seq.translate(str.maketrans("ATGCUN", "TACGAN"))
        return rev_seq
    elif file is not None:
        out_file = open("output_revcom.fasta", 'w')
        fasta_iter = fasta_reader(file)
        for record in fasta_iter:
            fasta_header, seq = record
            rev_seq = seq[::-1]
            rev_seq = rev_seq.translate(str.maketrans("ATGCUN", "TACGAN"))
            out_file.write(">" + fasta_header + "\n" + rev_seq + "\n")
        out_file.close()

# extract subseq from genome sequence
def ext_subseq(file="fasta_file", id="chr", st="start", end="end", strand="plus"):
    fasta_iter = fasta_reader(file)
    for record in fasta_iter:
        fasta_header, seq = record
        if id == fasta_header.strip() and strand == "plus":
            # -1 is necessary as it counts from 0
            sub_seq = seq[int(st-1):int(end)]
            print(sub_seq)
        elif id == fasta_header.strip() and strand == "minus":
            seq = rev_com(seq)
            sub_seq = seq[int(st-1):int(end)]
            print(sub_seq)

def fastq_format_check(file="fastq_file"):
    read_file = open(file, 'rU')
    x = 0
    for line in read_file:
        header = line.rstrip()
        if not header.startswith('@'):
            x = 1
        else:
            x = 0
        break
    return x

def tcsv(file="tab_file"):
    tab_file = csv.reader(open(file, 'r'), dialect=csv.excel_tab)
    csv_file = csv.writer(open('out.csv', 'w', newline=''), dialect=csv.excel)

    for record in tab_file:
        csv_file.writerow(record)

def ttsam(table="table", xfac=None, res=None, evar=True):
    d = pd.read_csv(table)
    if xfac and res is None:
        print("Error: xfac or res variable is missing")
        sys.exit(1)
    levels = d[xfac].unique()
    if len(levels) > 2:
        print("Error: there must be only two levels")
        sys.exit(1)
    a_val = d.loc[d[xfac] == levels[0], res].to_numpy()
    b_val = d.loc[d[xfac] == levels[1], res].to_numpy()
    a_count, b_count = len(a_val), len(b_val)
    count = [a_count, b_count]
    mean = d.groupby(xfac)[res].mean().to_numpy()
    sem = d.groupby(xfac)[res].sem().to_numpy()
    # degree of freedom
    # a_count, b_count = np.split(count, 2)
    dfa = a_count - 1
    dfb = b_count - 1
    # sample variance
    var_a = np.var(a_val, ddof=1)
    var_b = np.var(b_val, ddof=1)
    mean_diff = mean[0] - mean[1]
    # variable 95% CI
    varci_low = []
    varci_up = []
    tcritvar = [(stats.t.ppf((1 + 0.95) / 2, dfa)), (stats.t.ppf((1 + 0.95) / 2, dfb))]
    for i in range(len(levels)):
        varci_low.append(mean[i]-(tcritvar[i]*sem[i]))
        varci_up.append(mean[i]+(tcritvar[i]*sem[i]))

    var_test = 'equal'
    # perform levene to check for equal variances
    w, pvalue = stats.levene(a_val, b_val)
    if pvalue < 0.05:
        print(colored("Warning: the two group variance are not equal. Rerun the test with evar=False"))

    if evar is True:
        # pooled variance
        p_var = (dfa * var_a + dfb * var_b) / (dfa + dfb)
        # std error
        se = np.sqrt(p_var * (1.0 / a_count + 1.0 / b_count))
        df = dfa + dfb
    else:
        # Welch's t-test for unequal variance
        # calculate se
        a_temp = var_a/a_count
        b_temp = var_b/b_count
        df = ((a_temp + b_temp)**2) / ((a_temp**2)/(a_count-1)+(b_temp**2)/(b_count-1))
        se = np.sqrt(a_temp+b_temp)
        var_test = 'unequal'

    tval = np.divide(mean_diff, se)
    oneside_pval = stats.t.sf(np.abs(tval), df)
    twoside_pval = oneside_pval * 2
    # 95% CI for diff
    # 2.306 t critical at 0.05
    tcritdiff = stats.t.ppf((1 + 0.95) / 2, df)
    diffci_low = mean_diff-(tcritdiff*se)
    diffci_up = mean_diff+(tcritdiff*se)

    # print results
    print("\ntwo sample", levels, "t-test with", var_test, "variance", "\n")
    print(tabulate([["Mean diff", mean_diff], ["t", tval], ["std error", se], ["df", df],
                        ["P-value (one-tail)", oneside_pval], ["P-value (two-tail)", twoside_pval],
                        ["Lower 95%", diffci_low], ["Upper 95%", diffci_up]]), "\n")
    print("Parameter estimates\n")
    print(tabulate([[levels[0], count[0], mean[0], sem[0], varci_low[0], varci_low[1]], [levels[1], count[1],
                       mean[1], sem[1], varci_up[0], varci_up[1]]], headers=["Level", "Number", "Mean", "Std Error",
                       "Lower 95%", "Upper 95%"]), "\n")


    fig = plt.figure()
    d.boxplot(column=res, by=xfac, grid=False)
    plt.ylabel(res)
    plt.savefig('ttsam_boxplot.png', format='png', bbox_inches='tight', dpi=300)


def chisq(table="table"):
    d = pd.read_csv(table, index_col=0)
    tabulate_list = []
    chi_ps, p_ps, dof_ps, expctd_ps = stats.chi2_contingency(d.to_dict('split')['data'])
    tabulate_list.append(["Pearson", dof_ps, chi_ps, p_ps])
    chi_ll, p_ll, dof_ll, expctd_ll = stats.chi2_contingency(d.to_dict('split')['data'], lambda_="log-likelihood")
    tabulate_list.append(["Log-likelihood", dof_ll, chi_ll, p_ll])

    mosaic_dict = dict()
    m = d.to_dict('split')

    for i in range(d.shape[0]):
        for j in range(d.shape[1]):
            mosaic_dict[(m['index'][i], m['columns'][j])] = m['data'][i][j]

    print("\nChi-squared test\n")
    print(tabulate(tabulate_list, headers=["Test", "Df", "Chi-square", "P-value"]))
    print("\nExpected frequency counts\n")
    print(tabulate(expctd_ps, headers=d.to_dict('split')['columns'], showindex="always"))

    labels = lambda k: "" if mosaic_dict[k] != 0 else ""
    mosaic(mosaic_dict, labelizer=labels)
    plt.savefig('mosaic.png', format='png', bbox_inches='tight', dpi=300)


