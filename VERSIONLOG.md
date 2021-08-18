v2.0.6 has the following updates and changes (August 18, 2021)
- New function `Fasta.split_seq` added in `analys` module for splitting the sequence into smaller subsequences 
- 
v2.0.5 has the following updates and changes (August 16, 2021)
- New function `HtsAna.merge_featureCount` added in `analys` module for merging the counts for all samples
  obtained from featureCounts 

v2.0.4 has the following updates and changes (May 10, 2021)
- New class `Fasta` replaced with old `fasta`
- `Fasta.split_fasta` function added in `analys` module

v2.0.3 has the following updates and changes (April 14, 2021)
- New class `GeneExpression` and `General` added
- MA plot updated for adding the label names to genes and plot title  
- New theme option added to set new background styles to plots (now work only for MA plot)
- `gene_exp` class will be deprecated in future releases

v2.0.2 has the following updates and changes (March 13, 2021)
- PCA biplot updated for the y scale
- blood pressure dataset added

v2.0.1 has the following updates and changes (March 11, 2021)
- Dark theme added for image backgrounds

v2.0.0 has the following updates and changes (March 10, 2021)
- `analys.visuz.cluster.biplot` function updated to create the biplot where variables are more than observations

v1.0.9 has the following updates and changes (March 07, 2021)
- `analys.visuz.marker.mhat` function updated to handle the Fst values for Manhattan plot
- Boolean `log-scale` parameter added for choice of minus log10 conversion of <i>p</i> values
- New marker dataset with Fst values added. This dataset is provided by the Vincent Appiah, which is downloaded from 
  the Pf3K Project (pilot data release 5). This dataset can be accessed using `analys.get_data('fst').data`.

v1.0.8  has the following updates and changes (February 14, 2021)
- Function for regression metrics added (`bioinfokit.analys.stat.reg_metric`)
- It calculates Root Mean Square Error (RMSE), Mean squared error (MSE), Mean absolute error (MAE), 
  and Mean absolute percent error (MAPE) from regression fit

v1.0.7  has the following updates and changes (January 30, 2021)
- Plant species richness dataset added for regression analysis

v1.0.6  has the following updates and changes (January 29, 2021)
- Individual log fold change and p value cutoff added for volcano and inverted volcano plot 
  (`bioinfokit.visuz.gene_exp.volcano` and `bioinfokit.visuz.gene_exp.involcano`)

v1.0.5  has the following updates and changes (December 22, 2020)
- `analys.gff.gff_to_gtf` function updated to handle dot value for phase in CDS features

v1.0.4  has the following updates and changes (November 22, 2020)
- `Breast Cancer Wisconsin (Diagnostic) Data Set added
- `visuz.stat.roc` function added for visualizing the ROC

v1.0.3  has the following updates and changes (November 06, 2020)
- `bartlett` and `levene` function added to `analys.stat` class for checking the ANOVA assumptions
  for datasets in stacked format
- `tukey_hsd` function updated for grouping order
- Pandas series added as input for `fasta.extract_seq` function

v1.0.2  has the following updates and changes (October 26, 2020)
- `extract_seq` function moved to `fasta` class
- `extract_seq` function deprecated from `analys`
- visualization for single and multiple statistical bar charts updated for future releases
 
v1.0.1  has the following updates and changes (October 24, 2020)
- Tukey HSD test updated for interaction effect. Pairwise comparison for interaction effect can be calculated.
- `gff_to_gtf` function updated for the GFF3 file for non-coding RNA transcripts. GFF3 files with non-coding transcripts 
  (e.g. from miRBase GFF3) can be converted to GTF

v1.0.0  has the following updates and changes (October 10, 2020)
- genFam enrichment analysis function added (`bioinfokit.analys.genfam.fam_enrich`)
- genfam test added

v0.9.9  has the following updates and changes (October 04, 2020)
- Tukey HSD test added to perform multiple pairwise comparisons (`bioinfokit.analys.stat.tukey_hsd`)

v0.9.8  has the following updates and changes (September 25, 2020)
- new option `mrna_feature_name` added in `analys.gff.gff_to_gtf` if the name of the feature (column 3 of GFF3 file) of 
  protein coding mRNA is other than 'mRNA' or 'transcript' (e.g. some GFF3 file has this feature named as 
  protein_coding_gene )
-  `dim` option added to `visuz.cluster.screeplot`, `visuz.cluster.pcaplot` and `visuz.cluster.biplot` to control the
   figure size

v0.9.7  has the following updates and changes (September 18, 2020)
- `seqcov` moved to `fastq` class
-  `sra_db` function added under `fastq` class for batch download of FASTQ files 
  from NCBI SRA database 

v0.9.6  has the following updates and changes (August 22, 2020)
- In t-test, the one sample t and paired t-test added
- Two sample t-test switched to class based method
- t-test function name changed to `ttest` from `ttsam`
- programmatic access to chi-squared independence test dataset added
- boxplot removed from t-test
- 'adjustText' module added in `setup.py` (issue #12)

v0.9.5  has the following updates and changes (August 14, 2020)
- In chi-squared test, the sum of probabilities is rounded to 10 for exact sum in case of floats

v0.9.4  has the following updates and changes (August 13, 2020)
- chi-squared goodness of fit test added under the `stat.chisq`
- chi-squared independence test updated for output as class attributes and mosaic plot removed
- `mergevcf` renamed to `concatvcf` to keep with conventional naming (issue # 9)
- programmatic access to chi-squared independence test dataset added
- `marker.vcf_anot` function updated for tab-delimited text output

v0.9.3  has the following updates and changes (August 08, 2020)
- The error message for volcano, inverted volcano, and MA plot updated
  when there are no significant or non-significant genes (issue # 7)
- The `vcf_anot` function output updated for strand information 

v0.9.2  has the following updates and changes (July 30, 2020)
- The manhatten plot updated to add the lables in sorted order for numerical strings 
- The manhatten plot updated to add figname option 

v0.9.1  has the following updates and changes (July 30, 2020)
- TPM normalization function added

v0.9  has the following updates and changes (July 28, 2020)
- RPKM normalization function added

v0.8.9  has the following updates and changes (July 28, 2020)
- gene expression raw count normalization class added as 'analys.norm'
- CPM normalization function added

v0.8.8  has the following updates and changes (July 02, 2020)
- check for lfc_thr and pv_thr added

v0.8.7  has the following updates and changes (July 01, 2020)
- legend labels, position, and figname parameters added in volcano plot
- utility to check the non-numeric values added for `ma`, `volcano` and `involcano`
- plotlegend parameter added to `ma`

v0.8.6  has the following updates and changes (June 27, 2020)
- the parameter for log fold change threshold lines added in MA plot 
- legend labels, position, and figname parameters added in MA plot

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
