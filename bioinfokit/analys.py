from sklearn.decomposition import PCA
import pandas as pd
import re
import numpy as np
from bioinfokit.visuz import screeplot, pcaplot, general
from itertools import groupby, chain, combinations
import string
import sys, csv
import matplotlib.pyplot as plt
import scipy.stats as stats
from tabulate import tabulate
from statsmodels.graphics.mosaicplot import mosaic
from textwrap3 import wrap
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from statsmodels.formula.api import ols
import statsmodels.api as sm
from sklearn.linear_model import LinearRegression
from decimal import Decimal
from pathlib import Path
from sklearn.metrics import mean_squared_error
from collections import defaultdict


def seqcov(file="fastq_file", gs="genome_size"):
    x = fastq.fastq_format_check(file)
    if x == 1:
        print("Error: Sequences are not in fastq format")
        sys.exit(1)
    num_reads, total_len = fastq.fqreadcounter(file)
    # haploid genome_size must be in Mbp; convert in bp
    gs = gs * 1e6
    cov = round(float(total_len / gs), 2)
    print("Sequence coverage for", file, "is", cov)

def mergevcf(file="vcf_file_com_sep"):
    general.depr_mes("bioinfokit.analys.marker.mergevcf")

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
    general.depr_mes("bioinfokit.analys.fastq.fqreadcounter")


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
    general.depr_mes("bioinfokit.analys.fastq.fastq_format_check")


def tcsv(file="tab_file"):
    tab_file = csv.reader(open(file, 'r'), dialect=csv.excel_tab)
    csv_file = csv.writer(open('out.csv', 'w', newline=''), dialect=csv.excel)

    for record in tab_file:
        csv_file.writerow(record)


def ttsam(df='dataframe', xfac=None, res=None, evar=True):
    general.depr_mes("bioinfokit.visuz.stat.ttsam")


def chisq(table="table"):
    general.depr_mes("bioinfokit.visuz.stat.chisq")


class fasta:
    def __init__(self):
        pass

    # adapted from https://www.biostars.org/p/710/
    def fasta_reader(file="fasta_file"):
        read_file = open(file, "rU")
        fasta_iter = (rec[1] for rec in groupby(read_file, lambda line: line[0] == ">"))
        for record in fasta_iter:
            fasta_header = record.__next__()[1:].strip()
            fasta_header = re.split("\s+", fasta_header)[0]
            seq = "".join(s.strip() for s in fasta_iter.__next__())
            yield fasta_header, seq

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


class fastq:
    def __init__(self):
        pass

    def fastq_reader(file="fastq_file"):
        fastq_file = open(file, "r")
        for line in fastq_file:
            header_1 = line.rstrip()
            read = next(fastq_file).rstrip()
            header_2 = next(fastq_file).rstrip()
            read_qual_asc = next(fastq_file).rstrip()
            yield header_1, read, header_2, read_qual_asc

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
        num_reads = num_lines / 4
        return num_reads, total_len

    def fastq_format_check(file="fastq_file"):
        read_file = open(file, 'r')
        x = 0
        for line in read_file:
            header = line.rstrip()
            if not header.startswith('@'):
                x = 1
            else:
                x = 0
            break
        return x

    def detect_fastq_variant(file="fastq_file"):
        count = 0
        check = []
        fastq_file = open(file, 'rU')

        for line in fastq_file:
            header_1 = line.rstrip()
            read = next(fastq_file).rstrip()
            header_2 = next(fastq_file).rstrip()
            read_qual_asc = next(fastq_file).rstrip()
            asc_list = list(read_qual_asc)
            asc_list = list(map(ord, asc_list))
            min_q = min(asc_list)
            max_q = max(asc_list)
            check.append(min_q)
            check.append(max_q)
            count += 1
            if count == 40000:
                break
        fastq_file.close()
        min_q = min(check)
        max_q = max(check)
        if 64 > min_q >= 33 and max_q == 74:
            return 1
        elif min_q >= 64 and 74 < max_q <= 104:
            return 2
        elif 64 > min_q >= 33 and max_q <= 73:
            return 3

    def split_fastq(file="fastq_file"):
        x = fastq.fastq_format_check(file)
        if x == 1:
            print("Error: Sequences are not in sanger fastq format")
            sys.exit(1)
        fastq_iter = fastq.fastq_reader(file)
        out_file_name_1 = open(Path(file).stem+'_1.fastq', 'w')
        out_file_name_2 = open(Path(file).stem+'_2.fastq', 'w')
        i = 1
        for record in fastq_iter:
            header_1, read, header_2, read_qual_asc = record
            if (i % 2) == 0:
                out_file_name_2.write(header_1+'\n'+read+'\n'+header_2+'\n'+read_qual_asc+'\n')
            else:
                out_file_name_1.write(header_1+'\n'+read+'\n'+header_2+'\n'+read_qual_asc+'\n')
            i += 1

        out_file_name_1.close()
        out_file_name_2.close()


class marker:

    def __init__(self):
        pass

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

    def splitvcf(file='vcf_file', id='#CHROM'):
        read_vcf_file = open(file, 'r')
        info_lines, headers = [], []
        for line in read_vcf_file:
            if line.startswith(id):
                headers = line.strip().split('\t')
            elif line.startswith('##'):
                info_lines.append(line.strip())
        read_vcf_file.close()
        assert len(headers) != 0, "Non matching id parameter"
        read_vcf_file_df = pd.read_csv(file, sep='\t', comment='#', header=None)
        read_vcf_file_df.columns = headers
        chrom_ids = read_vcf_file_df[id].unique()
        for r in range(len(chrom_ids)):
            sub_df = read_vcf_file_df[read_vcf_file_df[id]==chrom_ids[r]]
            # out_vcf_file = open(chrom_ids[r]+'.vcf'
            with open(chrom_ids[r]+'.vcf', 'w') as out_vcf_file:
                for l in info_lines:
                    out_vcf_file.write(l+'\n')
            sub_df.to_csv(chrom_ids[r]+'.vcf', mode='a', sep='\t', index=False)
            out_vcf_file.close()

    def vcfreader(file='vcf_file', id='#CHROM'):
        read_vcf_file = open(file, 'r')
        info_lines, headers = [], []
        for line in read_vcf_file:
            if line.startswith(id):
                headers = line.strip().split('\t')
            elif line.startswith('##'):
                info_lines.append(line.strip())
            else:
                var_lines = line.strip().split('\t')
                yield headers, info_lines, var_lines
        read_vcf_file.close()
        assert len(headers) != 0, "Non matching id parameter"


    def vcf_anot(file='vcf_file', gff_file='gff_file', id='#CHROM', anot_attr=None):
        gff_iter = gff.gffreader(gff_file)
        gene_cord = defaultdict(list)
        cds_cord = defaultdict(list)
        exon_cord = defaultdict(list)
        ftr_cord = defaultdict(list)
        ttr_cord = defaultdict(list)
        sc_cord = defaultdict(list)
        st_cord = defaultdict(list)
        igenic_cord = defaultdict(list)
        intragenic_cord = defaultdict(list)
        # also for introns between the exons
        intragenic_cord_exon = defaultdict(list)
        gene_id_dict = dict()
        transcript_name_dict = dict()
        transcript_strand_dict = dict()
        chr_list = set([])
        for record in gff_iter:
            chr, gene_id, gene_name, transcript_id, source, feature_type, st, ende, strand, attr = record
            if feature_type == 'gene':
                if chr not in chr_list:
                    gene_number_1 = 1
                chr_list.add(chr)
                gene_cord[(chr, gene_id, gene_number_1)]=[st, ende]
                gene_id_dict[(chr, gene_number_1)] = gene_id
                gene_number_1 += 1
            elif feature_type == 'mRNA' or feature_type == 'transcript':
                cds_cord[(chr, transcript_id)] = []
                exon_cord[(chr, transcript_id)] = []
                ftr_cord[transcript_id] = []
                ttr_cord[transcript_id] = []
                sc_cord[transcript_id] = []
                st_cord[transcript_id] = []
                transcript_strand_dict[transcript_id] = strand
                if anot_attr:
                    transcript_name_dict[transcript_id] = re.search(anot_attr+'=(.+?)(;|$)',  attr).group(1)
            elif feature_type == 'CDS':
                cds_cord[(chr, transcript_id)].append([st, ende])
            elif feature_type == 'exon':
                exon_cord[(chr, transcript_id)].append([st, ende])
            elif feature_type == 'five_prime_UTR':
                ftr_cord[(chr, transcript_id)].append([st, ende])
            elif feature_type == 'three_prime_UTR':
                ttr_cord[(chr, transcript_id)].append([st, ende])
            elif feature_type == 'start_codon':
                sc_cord[(chr, transcript_id)].append([st, ende])
            elif feature_type == 'stop_codon':
                st_cord[(chr, transcript_id)].append([st, ende])

        # get intergenic regions
        for gene, cord in gene_cord.items():
            chr, gene_id, gene_number = gene[0], gene[1], gene[2]
            for x in chr_list:
                if x == chr and gene_number == 1:
                    igenic_cord[(chr, gene_id)] = [1, int(cord[0])-1]
                elif x == chr and gene_number != 1:
                    igenic_cord[(chr, gene_id)] = \
                        [int(gene_cord[(chr, gene_id_dict[(chr, int(gene_number)-1)], int(gene_number)-1)][1])+1, int(cord[0])-1]

        # get intragenic regions based on CDS
        for transcript, cord in cds_cord.items():
            chr, transcript_id = transcript[0], transcript[1]
            intragenic_cord[(chr, transcript_id)] = []
            for x in chr_list:
                if x == chr:
                    cord.sort(key=lambda k: k[0])
                    if len(cord) > 1:
                        for y in range(len(cord)-1):
                            intragenic_cord[(chr, transcript_id)].append([int(cord[y][1])+1, int(cord[y+1][0])-1])

        # get intragenic regions based on exon
        for transcript, cord in exon_cord.items():
            chr, transcript_id = transcript[0], transcript[1]
            intragenic_cord_exon[(chr, transcript_id)] = []
            for x in chr_list:
                if x == chr:
                    cord.sort(key=lambda k: k[0])
                    if len(cord) > 1:
                        for y in range(len(cord) - 1):
                            intragenic_cord_exon[(chr, transcript_id)].append([int(cord[y][1]) + 1, int(cord[y + 1][0]) - 1])

        def var_region_check(_dict, _chr,  _region, _anot_attr, _transcript_name_dict, _var_region, _transcript_name,
                             _transcript_id, _transcript_strand):
            for transcript, cord in _dict.items():
                for i in range(len(cord)):
                    if transcript[0] == chr and int(cord[i][0]) <= int(var_pos) <= int(cord[i][1]):
                        _var_region = _region
                        _transcript_id = transcript[1]
                        _transcript_strand = transcript_strand_dict[_transcript_id]
                        if anot_attr:
                            _transcript_name = transcript_name_dict[_transcript_id]
                        break
                if _var_region:
                    break
            return _var_region, _transcript_name, _transcript_id, _transcript_strand

        vcf_iter = marker.vcfreader(file, id)
        vcf_out_anot = open(Path(file).stem+'_anot.vcf', 'a')
        for_info_lines = 1
        transcript_id=None
        transcript_name=None
        transcript_strand=None
        transcript_name_return=transcript_name
        transcript_id_return=transcript_id
        transcript_strand_return=transcript_strand
        for record in vcf_iter:
            headers, info_lines, chr, var_pos = record[0], record[1], record[2][0], record[2][1]
            if for_info_lines == 1:
                for_info_lines = 0
                for l in info_lines:
                    vcf_out_anot.write(l+'\n')
                headers.extend(['genomic region', 'transcript ID', 'transcript name'])
                vcf_out_anot.write('\t'.join(x for x in headers) + '\n')

            var_region = None
            if var_region is None:
                for transcript, cord in igenic_cord.items():
                    if transcript[0] == chr and int(cord[0]) <= int(var_pos) <= int(cord[1]):
                        var_region = 'Intergenic'
                        transcript_id_return = None
                        transcript_strand_return = None
                        if anot_attr:
                            transcript_name_return = None
                        break
            if var_region is None:
                var_region, transcript_name_return, transcript_id_return, transcript_strand_return = var_region_check(cds_cord, chr, 'CDS',
                    anot_attr, transcript_name_dict, var_region, transcript_name, transcript_id, transcript_strand)
            if var_region is None:
                var_region, transcript_name_return, transcript_id_return, transcript_strand_return = var_region_check(ftr_cord, chr,
                    'five_prime_UTR', anot_attr, transcript_name_dict, var_region, transcript_name, transcript_id, transcript_strand)
            if var_region is None:
                var_region, transcript_name_return, transcript_id_return, transcript_strand_return = var_region_check(ttr_cord, chr,
                    'three_prime_UTR', anot_attr, transcript_name_dict, var_region, transcript_name, transcript_id, transcript_strand)
            if var_region is None:
                var_region, transcript_name_return, transcript_id_return, transcript_strand_return = var_region_check(sc_cord, chr,
                    'start_codon', anot_attr, transcript_name_dict, var_region, transcript_name, transcript_id, transcript_strand)
            if var_region is None:
                var_region, transcript_name_return, transcript_id_return, transcript_strand_return = var_region_check(st_cord, chr,
                    'stop_codon', anot_attr, transcript_name_dict, var_region, transcript_name, transcript_id, transcript_strand)
            if var_region is None:
                var_region, transcript_name_return, transcript_id_return, transcript_strand_return = var_region_check(exon_cord, chr,
                    'exon', anot_attr, transcript_name_dict, var_region, transcript_name, transcript_id, transcript_strand)

            '''
            if var_region is None:
                for transcript, cord in cds_cord.items():
                    for i in range(len(cord)):
                        if transcript[0] == chr and int(cord[i][0]) <= int(var_pos) <= int(cord[i][1]):
                            var_region = 'CDS'
                            transcript_id = transcript[1]
                            if anot_attr:
                                transcript_name = transcript_name_dict[transcript_id]
                            break
                    if var_region:
                        break
            
            if var_region is None:
                for transcript, cord in ftr_cord.items():
                    for i in range(len(cord)):
                        if transcript[0] == chr and int(cord[i][0]) <= int(var_pos) <= int(cord[i][1]):
                            var_region = 'five_prime_UTR'
                            break
                    if var_region:
                        break
                        
            if var_region is None:
                for transcript, cord in ttr_cord.items():
                    for i in range(len(cord)):
                        if transcript[0] == chr and int(cord[i][0]) <= int(var_pos) <= int(cord[i][1]):
                            var_region = 'three_prime_UTR'
                            break
                    if var_region:
                        break
                        
            if var_region is None:
                for transcript, cord in sc_cord.items():
                    for i in range(len(cord)):
                        if transcript[0] == chr and int(cord[i][0]) <= int(var_pos) <= int(cord[i][1]):
                            var_region = 'start_codon'
                            break
                    if var_region:
                        break
            if var_region is None:
                for transcript, cord in st_cord.items():
                    for i in range(len(cord)):
                        if transcript[0] == chr and int(cord[i][0]) <= int(var_pos) <= int(cord[i][1]):
                            var_region = 'stop_codon'
                            break
                    if var_region:
                        break            
            # keep exons at end as it also contains UTR part
            if var_region is None:
                for transcript, cord in exon_cord.items():
                    for i in range(len(cord)):
                        if transcript[0] == chr and int(cord[i][0]) <= int(var_pos) <= int(cord[i][1]):
                            var_region = 'exon'
                            break
                    if var_region:
                        break            
            '''
            if var_region is None:
                for transcript, cord in intragenic_cord.items():
                    transcript_strand_return = transcript_strand_dict[transcript[1]]
                    transcript_id_return = transcript[1]
                    if len(cord) >= 1:
                        for i in range(len(cord)):
                            if transcript[0] == chr and int(cord[i][0]) <= int(var_pos) <= int(cord[i][1]):
                                var_region = 'Introns'
                                break
                    if var_region:
                        break

            if var_region is None:
                for transcript, cord in intragenic_cord_exon.items():
                    transcript_strand_return = transcript_strand_dict[transcript[1]]
                    transcript_id_return = transcript[1]
                    if len(cord) >= 1:
                        for i in range(len(cord)):
                            if transcript[0] == chr and int(cord[i][0]) <= int(var_pos) <= int(cord[i][1]):
                                var_region = 'Introns'
                                break
                    if var_region:
                        break

            vcf_out_anot.write('\t'.join(str(x) for x in record[2])+'\t'+str(var_region)+'\t'+str(transcript_id_return)+
                               '\t'+str(transcript_name_return)+'\t'+str(transcript_strand_return)+'\n')


class format:
    def __init__(self):
        pass

    def fqtofa(file="fastq_file"):
        x = fastq.fastq_format_check(file)
        if x == 1:
            print("Error: Sequences are not in sanger fastq format")
            sys.exit(1)

        read_file = open(file, "rU")
        out_file = open("output.fasta", 'w')
        for line in read_file:
            header_1 = line.rstrip()
            read = next(read_file).rstrip()
            header_2 = next(read_file).rstrip()
            read_qual = next(read_file).rstrip()
            out_file.write(header_1+"\n"+'\n'.join(wrap(read, 60))+"\n")
        read_file.close()

    def tabtocsv(file="tab_file"):
        tab_file = csv.reader(open(file, 'r'), dialect=csv.excel_tab)
        csv_file = csv.writer(open('output.csv', 'w', newline=''), dialect=csv.excel)

        for record in tab_file:
            csv_file.writerow(record)

    def csvtotab(file="csv_file"):
        csv_file = csv.reader(open(file, 'r'), dialect=csv.excel)
        tab_file = csv.writer(open('output.txt', 'w', newline=''), dialect=csv.excel_tab)

        for record in csv_file:
            tab_file.writerow(record)

    def hmmtocsv(file="hmm_file"):
        hmm_file = open(file, "rU")
        csv_file = open("ouput_hmm.csv", "w")

        for line in hmm_file:
            line = line.strip()
            if not line.startswith("#"):
                data = re.split(' +', line)
                if len(data) == 19:
                    data[18] = data[18].replace(',', ' ')
                    csv_file.write(str.join(',', data))
                    csv_file.write("\n")
                elif len(data) > 19:
                    ele = list(range(18, len(data)))
                    data[18] = " ".join([e for i, e in enumerate(data) if i in ele])
                    data[18] = data[18].replace(',', '')
                    csv_file.write(str.join(',', data[0:19]))
                    csv_file.write("\n")
        hmm_file.close()
        csv_file.close()

    # find sanger fastq phred quality encoding format
    def fq_qual_var(file=None):
        if file is None:
            print("Error: No sanger fastq file provided")
            sys.exit(1)
        x = fastq.fastq_format_check(file)
        if x == 1:
            print("Error: Sequences are not in sanger fastq format")
            sys.exit(1)

        qual_format = fastq.detect_fastq_variant(file)

        if qual_format == 1:
            print("The fastq quality format is illumina 1.8+ (Offset +33)")
        elif qual_format == 2:
            print("The fastq quality format is illumina 1.3/1.4 (Offset +64)")
        elif qual_format == 3:
            print("The fastq quality format is Sanger (Offset +33)")
        else:
            print("\nError: Wrong quality format\n")
            sys.exit(1)


class stat:
    def __init__(self):
        pass

    def anova(self, df='dataframe', xfac=None, res=None):
        # drop NaN
        df = df.dropna()
        df = df[[xfac[0], res]]
        assert xfac and res is not None, "xfac or res variable is missing"
        grand_mean = df[res].mean()
        total_obs = df.count()[0]

        if len(xfac) == 1:
            levels = df[xfac[0]].unique()
            assert len(levels) > 2, 'levels must be more than 2; use two-sample t-test for two levels'
            levels.sort()
            ss_trt_between = np.sum(df.groupby(xfac).count() * (df.groupby(xfac).mean()-grand_mean)**2)[0]
            ss_err_within = 0
            for name, group in df.groupby(xfac):
                ss_err_within = ss_err_within + np.sum((group[res]-group[res].mean()) ** 2)
            ss_total = ss_trt_between + ss_err_within
            df_trt_between = len(levels)-1
            df_err_within = total_obs-len(levels)
            df_total = df_trt_between + df_err_within
            ms_trt_between = ss_trt_between / df_trt_between
            ms_err_within = ss_err_within / df_err_within
            f_value = ms_trt_between / ms_err_within
            p_value = '%.4E' % Decimal(stats.f.sf(f_value, df_trt_between, df_err_within))
            anova_table = []
            anova_table.append(
                ["Model", df_trt_between, ss_trt_between, round(ms_trt_between, 4), round(f_value, 4), p_value])
            anova_table.append(["Error", df_err_within, ss_err_within, round(ms_err_within, 4), "", ""])
            anova_table.append(["Total", df_total, ss_total, "", "", ""])
            print("\nANOVA Summary:\n")
            print(tabulate(anova_table, headers=["Source", "Df", "Sum Squares", "Mean Squares", "F", "Pr(>F)"]),
                  "\n")



    def oanova(table="table", res=None, xfac=None, ph=False, phalpha=0.05):
        # create and run model
        model = ols('{} ~ C({})'.format(res, xfac), data=table).fit()
        anova_table = sm.stats.anova_lm(model, typ=2)

        # treatments
        # this is for bartlett test
        levels = table[xfac].unique()
        fac_list = []
        data_summary = []
        for i in levels:
            temp_summary = []
            temp = table.loc[table[xfac]==i, res]
            fac_list.append(temp)
            temp_summary.append(i)
            temp_summary.extend(temp.describe().to_numpy())
            data_summary.append(temp_summary)

        print("\nTable Summary\n")
        print(tabulate(data_summary, headers=["Group", "Count", "Mean", "Std Dev", "Min", "25%", "50%", "75%", "Max"]), "\n")

        # check assumptions
        # Shapiro-Wilk  data is drawn from normal distribution.
        w, pvalue1 = stats.shapiro(model.resid)
        w, pvalue2 = stats.bartlett(*fac_list)
        if pvalue1 < 0.05:
            print("Warning: Data is not drawn from normal distribution")
        else:
            # samples from populations have equal variances.
            if pvalue2 < 0.05:
                print("Warning: treatments do not have equal variances")

        print("\nOne-way ANOVA Summary\n")
        print(anova_table)
        print("\n")

        # if post-hoc test is true
        if ph:
            # perform multiple pairwise comparison (Tukey HSD)
            m_comp = pairwise_tukeyhsd(endog=table[res], groups=table[xfac], alpha=phalpha)
            print("\nPost-hoc Tukey HSD test\n")
            print(m_comp, "\n")

        print("ANOVA Assumption tests\n")
        print("Shapiro-Wilk (P-value):", pvalue1, "\n")
        print("Bartlett (P-value):", pvalue2, "\n")

    def lin_reg(self, df="dataframe", y=None, x=None):
        df = df.dropna()
        assert x and y is not None, "Provide proper column names for X and Y variables"
        assert type(x) is list or type(y) is list, "X or Y column names should be list"
        # min data should be 4 or more
        assert df.shape[0] >= 4, "Very few data"
        self.X = df[x].to_numpy()
        self.Y = df[y].to_numpy()
        # number of independent variables
        p = len(x)
        # number of parameter estimates (+1 for intercept and slopes)
        e = p+1
        # number of samples/observations
        n = len(df[y])

        # run regression
        reg_out = LinearRegression().fit(self.X, self.Y)
        # coefficient  of determination
        r_sq = round(reg_out.score(self.X, self.Y), 4)
        # Correlation coefficient (r)
        # Adjusted r-Squared
        r_sq_adj = round(1 - (1 - r_sq) * ((n - 1)/(n-p-1)), 4)
        # RMSE
        # RMSE = standard deviation of the residuals
        rmse = round(np.sqrt(1-r_sq) * np.std(self.Y), 4)
        # intercept and slopes
        reg_intercept = reg_out.intercept_
        reg_slopes = reg_out.coef_
        # predicted values
        self.y_hat = reg_out.predict(self.X)
        # residuals
        self.residuals = self.Y - self.y_hat
        # sum of squares
        regSS = np.sum((self.y_hat - np.mean(self.Y)) ** 2)  # variation explained by linear model
        residual_sse = np.sum((self.Y - self.y_hat) ** 2)  # remaining variation
        sst = np.sum((self.Y - np.mean(self.Y)) ** 2)  # total variation

        eq = ""
        for i in range(p):
            eq = eq+' + '+ '(' + str(round(reg_slopes[0][i], 4))+'*'+x[i] + ')'

        self.reg_eq = str(round(reg_intercept[0], 4)) + eq

        # variance and std error
        # Residual variance = MSE and sqrt of MSE is res stnd error
        sigma_sq_hat = round(residual_sse/(n-e), 4)
        # residual std dev
        res_stdev = round(np.sqrt(sigma_sq_hat))
        # standardized residuals
        self.std_residuals = self.residuals/res_stdev

        # https://stackoverflow.com/questions/22381497/python-scikit-learn-linear-model-parameter-standard-error
        # std error
        X_mat = np.empty(shape=(n, e), dtype=np.float)
        X_mat[:, 0] = 1
        X_mat[:, 1:e] = self.X
        var_hat = np.linalg.inv(X_mat.T @ X_mat) * sigma_sq_hat
        standard_error = []
        for param in range(e):
            standard_error.append(round(np.sqrt(var_hat[param, param]), 4))

        # t = b1 / SE
        params = list(chain(*[["Intercept"], x]))
        estimates = list(chain(*[[reg_intercept[0]], reg_slopes[0]]))
        tabulate_list = []
        for param in range(e):
            tabulate_list.append([params[param], estimates[param], standard_error[param],
                                  estimates[param]/standard_error[param],
                                  '%.4E' % Decimal(stats.t.sf(np.abs(estimates[param]/standard_error[param]), n-e)*2)   ])

        # anova
        anova_table = []
        anova_table.append(["Model", p, regSS, round(regSS/p, 4), round((regSS/p)/(residual_sse/(n-e)), 4),
                            '%.4E' % Decimal(stats.f.sf((regSS/p)/(residual_sse/(n-e)), p, n-e))])
        anova_table.append(["Error", n-e, residual_sse, round(residual_sse/(n-e), 4), "", ""])
        anova_table.append(["Total", n-1, sst, "", "", ""])


        print("\nRegression equation:\n")
        print(self.reg_eq)
        print("\nRegression Summary:")
        print(tabulate([["Dependent variables", x], ["Independent variables", y],
                        ["Coefficient of determination (r-squared)", r_sq], ["Adjusted r-squared", r_sq_adj],
                        ["Root Mean Square Error (RMSE)", rmse],
                        ["Mean of Y", round(np.mean(self.Y), 4)], ["Residual standard error", round(np.sqrt(sigma_sq_hat), 4)],
                        ["No. of Observations", n]], "\n"))
        print("\nRegression Coefficients:\n")
        print(tabulate(tabulate_list, headers=["Parameter", "Estimate", "Std Error", "t-value", "P-value Pr(>|t|)"]), "\n")
        print("\nANOVA Summary:\n")
        print(tabulate(anova_table, headers=["Source", "Df", "Sum Squares", "Mean Squares", "F", "Pr(>F)"]),
              "\n")

        # VIF for MLR
        # VIF computed as regressing X on remaining X
        # using correlation
        if p > 1:
            vif_table = []
            vif_df = df[x]
            df_corr = vif_df.corr()
            vif_mat = np.linalg.inv(df_corr)
            self.vif = vif_mat.diagonal()
            for i in range(len(self.vif)):
                vif_table.append([x[i], self.vif[i]])
            print("\nVariance inflation factor (VIF)\n")
            print(tabulate(vif_table, headers=["Variable", "VIF"]),
                  "\n")

        '''
        vif = []
        for i in range(len(x)):
            temp = x[:]
            print(i, x, temp)
            vif_y = x[i]
            del temp[i]
            vif_x = temp
            y_mat = df[vif_y]
            x_mat = df[vif_x]
            print(y_mat, '\n', x_mat, '\n')
            vif_reg_out = LinearRegression().fit(x_mat, y_mat)
            vif.append(1 / (1-vif_reg_out.score(x_mat, y_mat)))

        print(vif)
        '''

    def ttsam(df='dataframe', xfac=None, res=None, evar=True, alpha=0.05):
        # drop NaN
        df = df.dropna()
        if xfac and res is None:
            raise Exception("xfac or res variable is missing")
        levels = df[xfac].unique()
        levels.sort()
        if len(levels) != 2:
            raise Exception("there must be only two levels")
        a_val = df.loc[df[xfac] == levels[0], res].to_numpy()
        b_val = df.loc[df[xfac] == levels[1], res].to_numpy()
        a_count, b_count = len(a_val), len(b_val)
        count = [a_count, b_count]
        mean = [df.loc[df[xfac] == levels[0], res].mean(), df.loc[df[xfac] == levels[1], res].mean()]
        sem = [df.loc[df[xfac] == levels[0], res].sem(), df.loc[df[xfac] == levels[1], res].sem()]
        sd = [df.loc[df[xfac] == levels[0], res].std(), df.loc[df[xfac] == levels[1], res].std()]
        ci = (1-alpha)*100
        # degree of freedom
        # a_count, b_count = np.split(count, 2)
        dfa = a_count - 1
        dfb = b_count - 1
        # sample variance
        with np.errstate(invalid='ignore'):
            var_a = np.nan_to_num(np.var(a_val, ddof=1))
            var_b = np.nan_to_num(np.var(b_val, ddof=1))
        mean_diff = mean[0] - mean[1]
        # variable 95% CI
        varci_low = []
        varci_up = []
        tcritvar = [(stats.t.ppf((1 + (1-alpha)) / 2, dfa)), (stats.t.ppf((1 + (1-alpha)) / 2, dfb))]
        for i in range(len(levels)):
            varci_low.append(mean[i] - (tcritvar[i] * sem[i]))
            varci_up.append(mean[i] + (tcritvar[i] * sem[i]))

        var_test = 'equal'
        # perform levene to check for equal variances
        w, pvalue = stats.levene(a_val, b_val)
        if pvalue < alpha:
            print("Warning: the two group variance are not equal. Rerun the test with evar=False")

        if evar is True:
            # pooled variance
            p_var = (dfa * var_a + dfb * var_b) / (dfa + dfb)
            # std error
            se = np.sqrt(p_var * (1.0 / a_count + 1.0 / b_count))
            dfr = dfa + dfb
        else:
            # Welch's t-test for unequal variance
            # calculate se
            if a_count == 1 or b_count == 1:
                raise Exception('Not enough observation for either levels. The observations should be > 1 for both levels')
            a_temp = var_a / a_count
            b_temp = var_b / b_count
            dfr = ((a_temp + b_temp) ** 2) / ((a_temp ** 2) / (a_count - 1) + (b_temp ** 2) / (b_count - 1))
            se = np.sqrt(a_temp + b_temp)
            var_test = 'unequal'

        tval = np.divide(mean_diff, se)
        oneside_pval = stats.t.sf(np.abs(tval), dfr)
        twoside_pval = oneside_pval * 2
        # 95% CI for diff
        # 2.306 t critical at 0.05
        tcritdiff = stats.t.ppf((1 + (1-alpha)) / 2, dfr)
        diffci_low = mean_diff - (tcritdiff * se)
        diffci_up = mean_diff + (tcritdiff * se)

        # print results
        print("\nTwo sample", levels, "t-test with", var_test, "variance", "\n")
        print(tabulate([["Mean diff", mean_diff], ["t", tval], ["Std Error", se], ["df", dfr],
                        ["P-value (one-tail)", oneside_pval], ["P-value (two-tail)", twoside_pval],
                        ["Lower "+str(ci)+"%", diffci_low], ["Upper "+str(ci)+"%", diffci_up]]), "\n")
        print("Parameter estimates\n")
        print(tabulate([[levels[0], count[0], mean[0], sd[0], sem[0], varci_low[0], varci_up[0]], [levels[1], count[1],
                                                                                             mean[1], sd[1], sem[1],
                                                                                             varci_low[1], varci_up[1]]],
                       headers=["Level", "Number", "Mean", "Std Dev", "Std Error",
                                "Lower "+str(ci)+"%", "Upper "+str(ci)+"%"]), "\n")

        fig = plt.figure()
        df.boxplot(column=res, by=xfac, grid=False)
        plt.ylabel(res)
        plt.savefig('ttsam_boxplot.png', format='png', bbox_inches='tight', dpi=300)

    def chisq(df='dataframe'):
        # d = pd.read_csv(table, index_col=0)
        tabulate_list = []
        chi_ps, p_ps, dof_ps, expctd_ps = stats.chi2_contingency(df.to_dict('split')['data'])
        tabulate_list.append(["Pearson", dof_ps, chi_ps, p_ps])
        chi_ll, p_ll, dof_ll, expctd_ll = stats.chi2_contingency(df.to_dict('split')['data'], lambda_="log-likelihood")
        tabulate_list.append(["Log-likelihood", dof_ll, chi_ll, p_ll])

        mosaic_dict = dict()
        m = df.to_dict('split')

        for i in range(df.shape[0]):
            for j in range(df.shape[1]):
                mosaic_dict[(m['index'][i], m['columns'][j])] = m['data'][i][j]

        print("\nChi-squared test\n")
        print(tabulate(tabulate_list, headers=["Test", "Df", "Chi-square", "P-value"]))
        print("\nExpected frequency counts\n")
        print(tabulate(expctd_ps, headers=df.to_dict('split')['columns'], showindex="always"))

        labels = lambda k: "" if mosaic_dict[k] != 0 else ""
        mosaic(mosaic_dict, labelizer=labels)
        plt.savefig('mosaic.png', format='png', bbox_inches='tight', dpi=300)


class gff:
    def __init__(self):
        pass

    def gff_to_gtf(file='gff_file'):
        read_gff_file_cds = open(file, 'r')
        cds_dict_st, cds_dict_st_phase  = dict(), dict()
        cds_dict_end, cds_dict_end_phase = dict(), dict()
        cds_ct = 0
        for line in read_gff_file_cds:
            if not line.startswith('#'):
                line = re.split('\s+', line.strip())
                if line[2]=='mRNA' or line[2]=='transcript':
                    # attr = re.split(';', line[8])
                    # transcript_id = attr[0].split('=')[1]
                    if 'ID=' in line[8]:
                        transcript_id = re.search('ID=(.+?)(;|$)', line[8]).group(1)
                    else:
                        raise Exception(
                                "ID required in GFF3 file in attribute field for mRNA/transcript"
                                " feature type")
                    cds_dict_st[transcript_id] = []
                    cds_dict_end[transcript_id] = []
                elif line[2] == 'CDS':
                    cds_ct += 1
                    cds_dict_st[transcript_id].append(line[3])
                    cds_dict_end[transcript_id].append(line[4])
                    cds_dict_st_phase[(transcript_id, line[3])] = line[7]
                    cds_dict_end_phase[(transcript_id, line[4])] = line[7]

        read_gff_file_cds.close()

        # check if CDS feature present in GFF3 file
        if cds_ct == 0:
            print ("Warning: No CDS feature type found in given GFF3 file. GTF file requires CDS feature type\n")
        read_gff_file = open(file, 'r')
        out_gtf_file = open(Path(file).stem+'.gtf', 'w')
        gene_id = ''
        transcript_id = ''
        gene_name = ''
        first_cds_present, last_cds_present, start_codon_present, end_codon_present = 0, 0, 0, 0

        for line in read_gff_file:
            if not line.startswith('#'):
                line = re.split('\s+', line.strip())
                if line[2]=='gene':
                    # attr = re.split(';', line[8])
                    if 'ID=' in line[8]:
                        gene_id = re.search('ID=(.+?)(;|$)',  line[8]).group(1)

                    if 'Name=' in line[8]:
                        gene_name = re.search('Name=(.+?)(;|$)',  line[8]).group(1)
                    elif 'gene_name=' in line[8]:
                        gene_name = re.search('gene_name=(.+?)(;|$)',  line[8]).group(1)
                    elif 'gene_id=' in line[8]:
                        gene_name = re.search('gene_id=(.+?)(;|$)',  line[8]).group(1)

                    if 'ID=' not in line[8]:
                        raise Exception("ID field required in GFF3 file in attribute field for gene feature type")

                    # gene_id = attr[0].split('=')[1]
                    # gene_attr_gtf = 'gene_id "'+gene_id+'"; gene_name "'+attr[1].split('=')[1]+'"; gene_source "'+\
                    #                line[1]+'";'
                    gene_attr_gtf = 'gene_id "' + gene_id + '"; gene_name "' + gene_name + '"; gene_source "' + line[1]+'";'
                    out_gtf_file.write('\t'.join(line[0:8])+'\t'+gene_attr_gtf+'\n')
                elif line[2]=='mRNA' or line[2]=='transcript':
                    cds_i, exon_i, ftr_i, ttr_i = 1, 1, 1, 1
                    # attr = re.split(';', line[8])
                    if 'ID=' in line[8]:
                        transcript_id = re.search('ID=(.+?)(;|$)', line[8]).group(1)
                    # if 'Parent=' in line[8]:
                    #    gene_id = re.search('Parent=(.*);', line[8]).group(1)

                    if 'ID=' not in line[8]:
                        raise Exception("ID field required in GFF3 file in attribute field for mRNA/transcript"
                                        " feature type")

                    # transcript_id = attr[0].split('=')[1]
                    # replace mRNA with transcript for gtf file
                    if line[2]=='mRNA':
                        line[2] = 'transcript'
                    # gene_attr_gtf = 'gene_id "' + gene_id + '"; transcript_id "' + transcript_id +'"; gene_name "'+\
                    #                attr[1].split('=')[1]+ '"; gene_source "' + line[1] + '";'
                    gene_attr_gtf = 'gene_id "' + gene_id + '"; transcript_id "' + transcript_id + \
                                    '"; gene_source "' + line[1] + '";'
                    out_gtf_file.write('\t'.join(line[0:8])+'\t'+gene_attr_gtf+'\n')
                elif line[2]=='CDS':
                    # attr = re.split(';', line[8])
                    if 'Parent=' in line[8]:
                        transcript_id_temp = re.search('Parent=(.+?)(;|$)', line[8]).group(1)

                    if 'Parent=' not in line[8]:
                        raise Exception(
                            "Parent field required in GFF3 file in attribute field for CDS"
                            " feature type")

                    # transcript_id_temp = attr[1].split('=')[1]

                    if line[3] == min(cds_dict_st[transcript_id_temp], key=int):
                        first_cds_present = 1
                    if line[4] == max(cds_dict_end[transcript_id_temp], key=int):
                        last_cds_present = 1

                    if transcript_id_temp == transcript_id:
                            # gene_attr_gtf = 'gene_id "' + gene_id + '"; transcript_id "' + transcript_id +'"; three_prime_UTR_number "'\
                            #                +str(ttr_i)+'"; gene_name "'+attr[1].split('=')[1]+ '"; gene_source "' + line[1] + '";'
                        gene_attr_gtf = 'gene_id "' + gene_id + '"; transcript_id "' + transcript_id + '"; cds_number "' \
                                        +str(cds_i)+'"; gene_name "'+gene_name+ '"; gene_source "' + line[1] + '";'
                        cds_i += 1
                    '''
                    if transcript_id_temp == transcript_id:
                        # gene_attr_gtf = 'gene_id "' + gene_id + '"; transcript_id "' + transcript_id +'"; cds_number "'\
                        #                +str(cds_i)+'"; gene_name "'+attr[1].split('=')[1]+ '"; gene_source "' + line[1] + '";'
                        gene_attr_gtf = 'gene_id "' + gene_id + '"; transcript_id "' + transcript_id + '"; cds_number "' \
                                            +str(cds_i)+'"; gene_source "' + line[1] + '";'
                        cds_i += 1
                    '''
                    out_gtf_file.write('\t'.join(line[0:8])+'\t'+gene_attr_gtf+'\n')
                elif line[2]=='exon':
                    # attr = re.split(';', line[8])
                    # transcript_id_temp = attr[1].split('=')[1]
                    if 'Parent=' in line[8]:
                        transcript_id_temp = re.search('Parent=(.+?)(;|$)', line[8]).group(1)

                    if  'Parent=' not in line[8]:
                        raise Exception(
                            "Parent field required in GFF3 file in attribute field for exon"
                            " feature type")

                    if transcript_id_temp == transcript_id:
                        # gene_attr_gtf = 'gene_id "' + gene_id + '"; transcript_id "' + transcript_id +'"; three_prime_UTR_number "'\
                        #                +str(ttr_i)+'"; gene_name "'+attr[1].split('=')[1]+ '"; gene_source "' + line[1] + '";'
                        gene_attr_gtf = 'gene_id "' + gene_id + '"; transcript_id "' + transcript_id + '"; exon_number "' \
                                        +str(exon_i)+'"; gene_name "'+gene_name+ '"; gene_source "' + line[1] + '";'
                        exon_i += 1

                    '''
                    if transcript_id_temp == transcript_id:
                        gene_attr_gtf = 'gene_id "' + gene_id + '"; transcript_id "' + transcript_id +'"; exon_number "'\
                                        +str(exon_i)+'"; gene_name "'+attr[1].split('=')[1]+ '"; gene_source "' + line[1] + '";'
                        exon_i += 1
                    '''
                    out_gtf_file.write('\t'.join(line[0:8])+'\t'+gene_attr_gtf+'\n')
                elif line[2]=='five_prime_UTR':
                    if 'Parent=' in line[8]:
                        transcript_id_temp = re.search('Parent=(.+?)(;|$)', line[8]).group(1)

                    if 'Parent=' not in line[8]:
                        raise Exception(
                            "Parent field required in GFF3 file in attribute field for five_prime_UTR"
                            " feature type")
                    # attr = re.split(';', line[8])
                    # transcript_id_temp = attr[1].split('=')[1]
                    if transcript_id_temp == transcript_id:
                        # gene_attr_gtf = 'gene_id "' + gene_id + '"; transcript_id "' + transcript_id +'"; three_prime_UTR_number "'\
                        #                +str(ttr_i)+'"; gene_name "'+attr[1].split('=')[1]+ '"; gene_source "' + line[1] + '";'
                        gene_attr_gtf = 'gene_id "' + gene_id + '"; transcript_id "' + transcript_id + '"; five_prime_UTR_number "' \
                                        +str(ftr_i)+'"; gene_name "'+gene_name+ '"; gene_source "' + line[1] + '";'
                        ftr_i += 1
                    '''
                    if transcript_id_temp == transcript_id:
                        gene_attr_gtf = 'gene_id "' + gene_id + '"; transcript_id "' + transcript_id +'"; five_prime_UTR_number "'\
                                        +str(ftr_i)+'"; gene_name "'+attr[1].split('=')[1]+ '"; gene_source "' + line[1] + '";'
                        ftr_i += 1
                    '''
                    out_gtf_file.write('\t'.join(line[0:8])+'\t'+gene_attr_gtf+'\n')
                elif line[2]=='three_prime_UTR':
                    if 'Parent=' in line[8]:
                        transcript_id_temp = re.search('Parent=(.+?)(;|$)', line[8]).group(1)

                    if 'Parent=' not in line[8]:
                        raise Exception(
                            "Parent field required in GFF3 file in attribute field for three_prime_UTR"
                            " feature type")

                    # attr = re.split(';', line[8])
                    # transcript_id_temp = attr[1].split('=')[1]
                    if transcript_id_temp == transcript_id:
                        # gene_attr_gtf = 'gene_id "' + gene_id + '"; transcript_id "' + transcript_id +'"; three_prime_UTR_number "'\
                        #                +str(ttr_i)+'"; gene_name "'+attr[1].split('=')[1]+ '"; gene_source "' + line[1] + '";'
                        gene_attr_gtf = 'gene_id "' + gene_id + '"; transcript_id "' + transcript_id + '"; three_prime_UTR_number "' \
                                        +str(ttr_i)+'"; gene_name "'+gene_name+ '"; gene_source "' + line[1] + '";'
                        ttr_i += 1
                    out_gtf_file.write('\t'.join(line[0:8])+'\t'+gene_attr_gtf+'\n')
                elif line[2]=='start_codon':
                    start_codon_present = 1
                    # attr = re.split(';', line[8])
                    # transcript_id_temp = attr[1].split('=')[1]
                    if 'Parent=' in line[8]:
                        transcript_id_temp = re.search('Parent=(.+?)(;|$)', line[8]).group(1)

                    if 'Parent=' not in line[8]:
                        raise Exception(
                            "Parent field required in GFF3 file in attribute field for three_prime_UTR"
                            " feature type")
                    if transcript_id_temp == transcript_id:
                        # gene_attr_gtf = 'gene_id "' + gene_id + '"; transcript_id "' + transcript_id + '"; gene_name "'+\
                        #                attr[1].split('=')[1]+ '"; gene_source "' + line[1] + '";'
                        gene_attr_gtf = 'gene_id "' + gene_id + '"; transcript_id "' + transcript_id + '"; gene_name "' + \
                                            gene_name+ '"; gene_source "' + line[1] + '";'
                    out_gtf_file.write('\t'.join(line[0:8])+'\t'+gene_attr_gtf+'\n')
                elif line[2]=='stop_codon':
                    end_codon_present = 1
                    # attr = re.split(';', line[8])
                    # transcript_id_temp = attr[1].split('=')[1]
                    if 'Parent=' in line[8]:
                        transcript_id_temp = re.search('Parent=(.+?)(;|$)', line[8]).group(1)

                    if 'Parent=' not in line[8]:
                        raise Exception(
                            "Parent field required in GFF3 file in attribute field for three_prime_UTR"
                            " feature type")
                    if transcript_id_temp == transcript_id:
                        # gene_attr_gtf = 'gene_id "' + gene_id + '"; transcript_id "' + transcript_id + '"; gene_name "'+\
                        #                attr[1].split('=')[1]+ '"; gene_source "' + line[1] + '";'
                        gene_attr_gtf = 'gene_id "' + gene_id + '"; transcript_id "' + transcript_id + '"; gene_name "' + \
                                            gene_name+ '"; gene_source "' + line[1] + '";'
                    out_gtf_file.write('\t'.join(line[0:8])+'\t'+gene_attr_gtf+'\n')

                if first_cds_present == 1 and start_codon_present == 0:
                    first_cds_present = 0
                    if 'Parent=' in line[8]:
                        gene_id_temp = re.search('Parent=(.+?)(;|$)', line[8]).group(1)

                    if 'Parent=' not in line[8]:
                        raise Exception(
                            "Parent field required in GFF3 file in attribute field for CDS"
                            " feature type")
                    # attr = re.split(';', line[8])
                    # transcript_id_temp = attr[1].split('=')[1]
                    if transcript_id_temp == transcript_id or gene_id_temp == gene_id:
                        if line[6] == '+':
                            codon_min_cord = int(min(cds_dict_st[transcript_id_temp], key=int))
                            cds_phase = int(
                                cds_dict_st_phase[(transcript_id_temp, min(cds_dict_st[transcript_id_temp], key=int))])
                            line[2], line[3], line[4] = 'start_codon', codon_min_cord + cds_phase, \
                                                        codon_min_cord + cds_phase + 2
                        elif line[6] == '-':
                            codon_max_cord = int(max(cds_dict_end[transcript_id_temp], key=int))
                            cds_phase = int(
                                cds_dict_end_phase[(transcript_id_temp, max(cds_dict_end[transcript_id_temp], key=int))])
                            line[2], line[3], line[4] = 'start_codon', codon_max_cord - 2 - cds_phase, \
                                                        codon_max_cord - cds_phase
                        # gene_attr_gtf = 'gene_id "' + gene_id + '"; transcript_id "' + transcript_id + '"; gene_name "' + \
                        #                attr[1].split('=')[1] + '"; gene_source "' + line[1] + '";'
                        gene_attr_gtf = 'gene_id "' + gene_id + '"; transcript_id "' + transcript_id + '"; gene_name "' + \
                                            gene_name + '"; gene_source "' + line[1] + '";'
                        out_gtf_file.write('\t'.join(str(x) for x in line[0:8]) + '\t' + gene_attr_gtf + '\n')

                if last_cds_present == 1 and end_codon_present == 0:
                    last_cds_present = 0
                    if 'Parent=' in line[8]:
                        gene_id_temp = re.search('Parent=(.+?)(;|$)', line[8]).group(1)

                    if 'Parent=' not in line[8]:
                        raise Exception(
                            "Parent field required in GFF3 file in attribute field for CDS"
                            " feature type")

                    # attr = re.split(';', line[8])
                    # transcript_id_temp = attr[1].split('=')[1]
                    if transcript_id_temp == transcript_id or gene_id_temp == gene_id:
                        if line[6] == '+':
                            codon_max_cord = int(max(cds_dict_end[transcript_id_temp], key=int))
                            cds_phase = int(cds_dict_end_phase[(transcript_id_temp, max(cds_dict_end[transcript_id_temp], key=int))])
                            line[2], line[3], line[4] = 'stop_codon', codon_max_cord - 2, codon_max_cord
                        elif line[6] == '-':
                            codon_min_cord = int(min(cds_dict_st[transcript_id_temp], key=int))
                            cds_phase = int(cds_dict_st_phase[(transcript_id_temp, max(cds_dict_st[transcript_id_temp], key=int))])
                            line[2], line[3], line[4] = 'stop_codon', codon_min_cord, codon_min_cord + 2
                        # gene_attr_gtf = 'gene_id "' + gene_id + '"; transcript_id "' + transcript_id + '"; gene_name "' + \
                        #                attr[1].split('=')[1] + '"; gene_source "' + line[1] + '";'
                        gene_attr_gtf = 'gene_id "' + gene_id + '"; transcript_id "' + transcript_id + '"; gene_name "' + \
                                            gene_name + '"; gene_source "' + line[1] + '";'
                        out_gtf_file.write('\t'.join(str(x) for x in line[0:8]) + '\t' + gene_attr_gtf + '\n')

        read_gff_file.close()
        out_gtf_file.close()

    def gffreader(file='gff_file'):
        read_gff_file = open(file, 'r')
        transcript_id = ''
        for line in read_gff_file:
            if not line.startswith('#'):
                line = re.split('\t', line.strip())
                if line[2]=='gene':
                    if 'ID=' in line[8]:
                        gene_id = re.search('ID=(.+?)(;|$)',  line[8]).group(1)

                    if 'Name=' in line[8]:
                        gene_name = re.search('Name=(.+?)(;|$)',  line[8]).group(1)
                    elif 'gene_name=' in line[8]:
                        gene_name = re.search('gene_name=(.+?)(;|$)',  line[8]).group(1)
                    elif 'gene_id=' in line[8]:
                        gene_name = re.search('gene_id=(.+?)(;|$)',  line[8]).group(1)

                    if 'ID=' not in line[8]:
                        raise Exception("ID field required in GFF3 file in attribute field for gene feature type")
                    yield (line[0], gene_id, gene_name, transcript_id, line[1], line[2], line[3], line[4], line[6], line[8])
                elif line[2]=='mRNA' or line[2]=='transcript':
                    if 'ID=' in line[8]:
                        transcript_id = re.search('ID=(.+?)(;|$)', line[8]).group(1)
                    else:
                        raise Exception("ID field required in GFF3 file in attribute field for mRNA/transcript"
                                        " feature type")
                    yield (line[0], gene_id, gene_name, transcript_id, line[1], line[2], line[3], line[4], line[6], line[8])
                elif line[2]=='CDS':
                    if 'Parent=' in line[8]:
                        transcript_id_temp = re.search('Parent=(.+?)(;|$)', line[8]).group(1)
                    else:
                        raise Exception(
                            "Parent field required in GFF3 file in attribute field for CDS"
                            " feature type")
                    if transcript_id_temp == transcript_id:
                        yield (line[0], gene_id, gene_name, transcript_id, line[1], line[2], line[3], line[4], line[6], line[8])
                elif line[2]=='exon':
                    if 'Parent=' in line[8]:
                        transcript_id_temp = re.search('Parent=(.+?)(;|$)', line[8]).group(1)
                    else:
                        raise Exception(
                            "Parent field required in GFF3 file in attribute field for exon"
                            " feature type")
                    if transcript_id_temp == transcript_id:
                        yield (line[0], gene_id, gene_name, transcript_id, line[1], line[2], line[3], line[4], line[6], line[8])
                elif line[2]=='five_prime_UTR':
                    if 'Parent=' in line[8]:
                        transcript_id_temp = re.search('Parent=(.+?)(;|$)', line[8]).group(1)
                    else:
                        raise Exception(
                            "Parent field required in GFF3 file in attribute field for five_prime_UTR"
                            " feature type")
                    if transcript_id_temp == transcript_id:
                        yield (line[0], gene_id, gene_name, transcript_id, line[1], line[2], line[3], line[4], line[6], line[8])
                elif line[2]=='three_prime_UTR':
                    if 'Parent=' in line[8]:
                        transcript_id_temp = re.search('Parent=(.+?)(;|$)', line[8]).group(1)
                    else:
                        raise Exception(
                            "Parent field required in GFF3 file in attribute field for three_prime_UTR"
                            " feature type")
                    if transcript_id_temp == transcript_id:
                        yield (line[0], gene_id, gene_name, transcript_id, line[1], line[2], line[3], line[4], line[6], line[8])
                elif line[2]=='start_codon':
                    if 'Parent=' in line[8]:
                        transcript_id_temp = re.search('Parent=(.+?)(;|$)', line[8]).group(1)
                    else:
                        raise Exception(
                            "Parent field required in GFF3 file in attribute field for three_prime_UTR"
                            " feature type")
                    if transcript_id_temp == transcript_id:
                        yield (line[0], gene_id, gene_name, transcript_id, line[1], line[2], line[3], line[4], line[6], line[8])

                elif line[2]=='stop_codon':
                    if 'Parent=' in line[8]:
                        transcript_id_temp = re.search('Parent=(.+?)(;|$)', line[8]).group(1)
                    else:
                        raise Exception(
                            "Parent field required in GFF3 file in attribute field for three_prime_UTR"
                            " feature type")
                    if transcript_id_temp == transcript_id:
                        yield (line[0], gene_id, gene_name, transcript_id, line[1], line[2], line[3], line[4], line[6], line[8])
        read_gff_file.close()

class get_data:
    def __init__(self, data=None):
        if data=='mlr':
            self.data = pd.read_csv("https://reneshbedre.github.io/assets/posts/reg/test_reg.csv")
        elif data=='boston':
            self.data = pd.read_csv("https://reneshbedre.github.io/assets/posts/reg/boston.csv")
        elif data=='volcano':
            self.data = pd.read_csv("https://reneshbedre.github.io/assets/posts/volcano/testvolcano.csv")
        elif data=='ma':
            self.data = pd.read_csv("https://reneshbedre.github.io/assets/posts/ma/ma.csv")
        elif data=='hmap':
            self.data = pd.read_csv("https://reneshbedre.github.io/assets/posts/heatmap/hm_cot.csv")
        elif data=='mhat':
            self.data = pd.read_csv("https://reneshbedre.github.io/assets/posts/mhat/gwas_res_sim.csv")
        elif data=='bdot':
            self.data = pd.read_csv("https://reneshbedre.github.io/assets/posts/bardot/bardot.txt", sep="\t")
        elif data=='corr':
            self.data = pd.read_csv("https://reneshbedre.github.io/assets/posts/corr/corr_dataset.csv")
        elif data=='slr':
            self.data = pd.read_csv("https://reneshbedre.github.io/assets/posts/reg/test_reg_uni.csv")
        elif data=='ttest':
            self.data = pd.read_csv("https://reneshbedre.github.io/assets/posts/ttest/genotype.csv")
        elif data=='gexp':
            self.data = pd.read_csv("https://reneshbedre.github.io/assets/posts/pca/cot_pca.csv")
        elif data=='iris':
            self.data = pd.read_csv("https://reneshbedre.github.io/assets/posts/pca/iris.csv")
        elif data=='digits':
            self.data = pd.read_csv("https://reneshbedre.github.io/assets/posts/tsne/digits.csv")
        elif data=='pbmc':
            self.data = pd.read_csv("https://reneshbedre.github.io/assets/posts/tsne/pbmc_seurat_processes.csv")
        elif data=='ath_root':
            self.data = pd.read_csv("https://reneshbedre.github.io/assets/posts/tsne/ath_root_sub_seurat_processes.csv")
        else:
            print("Error: Provide correct parameter for data\n")


