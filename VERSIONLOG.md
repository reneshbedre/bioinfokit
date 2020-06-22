v0.8.5  has the following updates and changes (June 22, 2020)
- `tsneplot` added for t-SNE visualization
- in `bardot` drop NA value function added to ignore missing values to plot dots
- scRNA-seq dataset added (PBMC and Arabidopsis root cells)

v0.8.4  has the following updates and changes (June 17, 2020)
- `fasta_reader` and `rev_com` moved to newly created `fasta` class
- `tsneplot` and  `vcf_anot` initialized for future release
- more parameters added in `biplot` (cluster coloring, datapoints)
- `figname` added in `hmap` 

v0.8.3  has the following updates and changes (June 03, 2020)
- `ma` function updated for absolute expression counts
- `svg` figures added
- `pca` function will be deprecated in future release

v0.8.1 and v0.8.2  has the following updates and changes (May 31, 2020)
- 2D and 3D loadings plot, biplot and scree plot functions added under the
  `cluster` class for PCA
- programmatic access to iris and cotton dataset added
- `pca` function will be deprecated in future release

v0.8 has the following updates and changes (May 24, 2020)
- GFF3 to GTF file conversion utility added and updated under class `gff`

v0.7.3 has the following updates and changes (May 14, 2020)
- In manhatten plot (`visuz.marker.mhat`), the labeling issue with `markernames` parameter corrected (see issue # 4 on github for details;
  thanks to mkirchler for reporting)
- `gstyle` parameter added in manhatten plot for box style annotation

v0.7.2 has the following updates and changes (May 08, 2020)
- `splitvcf` function added for splitting VCF file into individual VCF files for each chromosome
- `mergevcf` moved to `analys.marker` class

v0.7.1 has the following updates and changes (April 24, 2020)
- `reg_lin` function updated for multiple regression
- degree of freedom fixed for t-test for regression coefficients
- VIF calculation for MLR updated
- functions `fastq_reader` and `fqreadcounter` moved to `fastq` class

v0.7 has the following updates and changes
- `split_fastq` function added for splitting individual (left and right) paired-end fastq files
  from single interleaved paired-end file
- GFF3 to GTF file conversion utility added under class `gff`
- two-sample and Welch's t-test updated for CI and alpha parameter added
- module termcolor removed
- Programmatic access of dataset for `ttsam` added

v0.6 has the following updates and changes
- Programmatic access of dataset added (class `get_data`)
- More features for figures added (`figtype`, `axtickfontsize`, `axtickfontname`, `axxlabel`, `axylabel`, `xlm`, `ylm`,
  `yerrlw`, `yerrcw`)
- In volcano plot, the typo for xlabel corrected (-log2(FoldChange) to log2(FoldChange))
- `help` will be deprecated in future release
- VIF calculation for MLR updated
- adjustText removed

v0.5 has the following updates and changes
- Linear regression analysis added in `analys.stat` class
- `volcano`, `involcano`, `ma` and `heatmap` functions moved to new `visuz.gen_exp` class
- In `volcano`, parameters for new box type labelling and threshold grid lines added
- `corr_mat` updated for new colormaps and moved to stat class
- To visualize the graph in console itself (e.g. Jupyter notebook), show parameter added
- Pandas dataframe input added for `volcano`, `involcano`, `corr_mat`, `ma`, `ttsam`, and `chisq`
- `ttsam` and `chisq` moved to `analys.stat` class
- graph control parameters added for  `volcano`, `involcano`, `ma`, and `heatmap`
- documentation can also be accessed at https://reneshbedre.github.io/blog/howtoinstall.html
- `help` will be deprecated in future release
- fixed the numpy bug in `visuz.stat.bardot`. The `int` cast added to generate number of samples, which does not accept
  float (See details of numpy bug: https://github.com/numpy/numpy/issues/15345)

v0.4 has the following updates and changes
- function `analyis.format.fq_qual_var()` added for detecting the FASTQ quality encoding format
- `help` module added command-line help message
- class `fastq` added for FASTQ related functions


v0.3 has the following updates and changes
- bar-dot plot function added
- command-line help message class added
- [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3698146.svg)](https://doi.org/10.5281/zenodo.3698146)
