from sklearn.decomposition import PCA
import pandas as pd
import numpy as np
from bioinfokit.visuz import screeplot, pcaplot


def seqcov(file="fastq_file", gs="genome_size"):
    num_reads, total_len = fqreadcounter(file)
    # haploid genome_size must be in Mbp; convert in bp
    gs = gs * 1e6
    cov = round(float(total_len / gs), 4)
    print(file, cov)


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