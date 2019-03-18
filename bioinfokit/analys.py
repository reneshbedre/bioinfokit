

def seqcov(file="fastq_file", gs="genome_size"):
    num_reads, total_len = fqreadcounter(file)
    # haploid genome_size must be in Mbp; convert in bp
    gs = gs * 1e-6
    cov = float(total_len/gs)
    print(file, cov, "X")


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


def fqreadcounter(file="fastq_file"):
    read_file = open(file, "rU")
    num_lines = 0
    total_len = 0
    for line in read_file:
        num_lines =+ 1
        header_1 = line.rstrip()
        read = next(read_file).rstrip()
        len_read = len(read)
        total_len =+ len_read
        header_2 = next(read_file).rstrip()
        read_qual = next(read_file).rstrip()
    read_file.close()
    num_reads = num_lines/4
    return num_reads, total_len