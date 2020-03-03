Bioinformatics data analysis and visualization toolkit


**<span style="color:#33a8ff">How to install:</span>**

bioinfokit requires
- Python 3
- NumPy 
- scikit-learn
- seaborn
- pandas
- matplotlib
- SciPy 

```
git clone https://github.com/reneshbedre/bioinfokit.git
cd bioinfokit
python3 setup.py install
```

<!--
From Python package index,

`pip install -i https://test.pypi.org/simple/ bioinfokit `

**<span style="color:#33a8ff">Functions:</span>**
-->


<b>Volcano plot</b>

`bioinfokit.visuz.volcano(table, lfc, pv, lfc_thr, pv_thr, color, valpha, geneid, genenames, gfont)`

Parameters | Description
------------ | -------------
`table` |Comma separated (csv) gene expression table having atleast gene IDs, log fold change, P-values or adjusted P-values columns
`lfc` | Name of a column having log fold change values [string][default:logFC]
`pv` | Name of a column having P-values or adjusted P-values [string][default:p_values]
`lfc_thr` | Log fold change cutoff for up and downregulated genes [float][default:1.0]
`pv_thr` | P-values or adjusted P-values cutoff for up and downregulated genes [float][default:0.05]
`color` | Tuple of two colors [tuple][default: ("green", "red")]
`valpha` | Transparency of points on volcano plot [float (between 0 and 1)][default: 1.0]
`geneid` | Name of a column having gene Ids. This is necessary for plotting gene label on the points [string][default: None]
`genenames` | Tuple of gene Ids to label the points. The gene Ids must be present in the geneid column. If this option set to "deg" it will label all genes defined by lfc_thr and pv_thr [string, tuple, dict][default: None]
`gfont` | Font size for genenames [float][default: 10.0]

Returns:

Volcano plot image in same directory (volcano.png)

<a href="https://reneshbedre.github.io/blog/volcano.html" target="_blank">Working example</a>



<b>MA plot</b>

`bioinfokit.visuz.ma(table, lfc, ct_count, st_count, pv_thr)`

Parameters | Description
------------ | -------------
`table` |Comma separated (csv) gene expression table having atleast gene IDs, log fold change, and counts (control and treatment) columns
`lfc` | Name of a column having log fold change values [default:logFC]
`ct_count` | Name of a column having count values for control sample [default:value1]
`st_count` | Name of a column having count values for treatment sample [default:value2]
`lfc_thr` | Log fold change cutoff for up and downregulated genes [default:1]


Returns:

MA plot image in same directory (ma.png)

<a href="https://reneshbedre.github.io/blog/ma.html" target="_blank">Working example</a>

<b>Inverted Volcano plot</b>

`bioinfokit.visuz.involcano(table, lfc, pv, lfc_thr, pv_thr, color, valpha, geneid, genenames, gfont)`

Parameters | Description
------------ | -------------
`table` |Comma separated (csv) gene expression table having atleast gene IDs, log fold change, P-values or adjusted P-values
`lfc` | Name of a column having log fold change values [default:logFC]
`pv` | Name of a column having P-values or adjusted P-values [default:p_values]
`lfc_thr` | Log fold change cutoff for up and downregulated genes [default:1]
`pv_thr` | P-values or adjusted P-values cutoff for up and downregulated genes [default:0.05]
`color` | Tuple of two colors [tuple][default: ("green", "red")]
`valpha` | Transparency of points on volcano plot [float (between 0 and 1)][default: 1.0]
`geneid` | Name of a column having gene Ids. This is necessary for plotting gene label on the points [string][default: None]
`genenames` | Tuple of gene Ids to label the points. The gene Ids must be present in the geneid column. If this option set to "deg" it will label all genes defined by lfc_thr and pv_thr [string, tuple, dict][default: None]
`gfont` | Font size for genenames [float][default: 10.0]

Returns:

Inverted volcano plot image in same directory (involcano.png)

<a href="https://reneshbedre.github.io/blog/volcano.html" target="_blank">Working example</a>


<b>Correlation matrix plot</b>

`bioinfokit.visuz.corr_mat(table, corm)`

Parameters | Description
------------ | -------------
`table` | Dataframe object with numerical variables (columns) to find correlation. Ideally, you should have three or more variables. Dataframe should not have identifier column.
`corm` | Correlation method [pearson,kendall,spearman] [default:pearson]


Returns:

Correlation matrix plot image in same directory (corr_mat.png)

<a href="https://reneshbedre.github.io/blog/corr.html" target="_blank">Working example</a>

<b>Merge VCF files</b>

`bioinfokit.analys.mergevcf(file)`

Parameters | Description
------------ | -------------
`file` | Multiple vcf files and separate them by comma

Returns:

Merged VCF file (merge_vcf.vcf)

<a href="https://reneshbedre.github.io/blog/mergevcf.html" target="_blank">Working example</a>


<b>Merge VCF files</b>

`bioinfokit.analys.mergevcf(file)`

Parameters | Description
------------ | -------------
`file` | Multiple vcf files and separate them by comma

Returns:

Merged VCF file (merge_vcf.vcf)

<a href="https://reneshbedre.github.io/blog/mergevcf.html" target="_blank">Working example</a>


<b>PCA</b>

`bioinfokit.analys.pca(table)`

Parameters | Description
------------ | -------------
`table` | Dataframe object with numerical variables (columns). Dataframe should not have identifier column.

Returns:

PCA summary, scree plot (screepca.png), and 2D/3D pca plots (pcaplot_2d.png and pcaplot_3d.png)

<a href="https://reneshbedre.github.io/blog/pca_3d.html" target="_blank">Working example</a>


<b>Reverse complement of DNA sequence</b>

`bioinfokit.analys.rev_com(sequence)`

Parameters | Description
------------ | -------------
`seq` | DNA sequence to perform reverse complement
`file` | DNA sequence in a fasta file

Returns:

Reverse complement of original DNA sequence

<a href="https://reneshbedre.github.io/blog/revcom.html" target="_blank">Working example</a>

<b>Sequencing coverage</b>

`bioinfokit.analys.seqcov(file, gs)`

Parameters | Description
------------ | -------------
`file` | FASTQ file
`gs` | Genome size in Mbp

Returns:

Sequencing coverage of the given FASTQ file

<a href="https://reneshbedre.github.io/blog/seqcov.html" target="_blank">Working example</a>

<b>Convert TAB to CSV file</b>

`bioinfokit.analys.tcsv(file)`

Parameters | Description
------------ | -------------
`file` | TAB delimited text file

Returns:

CSV delimited file (out.csv)

<b>Heatmap</b>

`bioinfokit.visuz.hmap(table, cmap='seismic', scale=True, dim=(6, 8), 
    clus=True, zscore=None, xlabel=True, ylabel=True, tickfont=(12, 12))`

Parameters | Description
------------ | -------------
`file` | CSV delimited data file. It should not have NA or missing values
`cmap` | Color Palette for heatmap [string][default: 'seismic']
`scale` | Draw a color key with heatmap [boolean (True or False)][default: True]
`dim` | heatmap figure size [tuple of two floats (width, height) in inches][default: (6, 8)]
`clus` | Draw hierarchical clustering with heatmap [boolean (True or False)][default: True]
`zscore` | Z-score standardization of row (0) or column (1). It works when clus is True. [None, 0, 1][default: None]
`xlable` | Plot X-label [boolean (True or False)][default: True]
`ylable` | Plot Y-label [boolean (True or False)][default: True]
`tickfont` | Fontsize for X and Y-axis tick labels [tuple of two floats][default: (14, 14)]

Returns:

heatmap plot (heatmap.png, heatmap_clus.png)

<a href="https://reneshbedre.github.io/blog/hmap.html" target="_blank">Working example</a>

<b>Venn Diagram</b>

`bioinfokit.visuz.venn(vennset, venncolor, vennalpha, vennlabel)`

Parameters | Description
------------ | -------------
`vennset` | Venn dataset for 3 and 2-way venn. Data should be in the format of (100,010,110,001,101,011,111) for 3-way venn and 2-way venn (10, 01, 11) [default: (1,1,1,1,1,1,1)]
`venncolor` | Color Palette for Venn [color code][default: ('#00909e', '#f67280', '#ff971d')]
`vennalpha` | Transparency of Venn  [float (0 to 1)][default: 0.5]
`vennlabel` | Labels to Venn [string][default: ('A', 'B', 'C')]

Returns:

Venn plot (venn3.png, venn2.png)

<a href="https://reneshbedre.github.io/blog/venn.html" target="_blank">Working example</a>

<b>Two sample t-test with equal and unequal variance</b>

`bioinfokit.analys.ttsam(table, xfac, res, evar)`

Parameters | Description
------------ | -------------
`table` | CSV delimited data file. It should be stacked table with independent (xfac) and dependent (res) variable columns.
`xfac` | Independent group column name with two levels [string][default: None]
`res` | Response variable column name [string][default: None]
`evar` | t-test with equal variance [bool (True or False)][default: True]

Returns:

summary output and group boxplot (ttsam_boxplot.png)

<a href="https://reneshbedre.github.io/blog/ttest.html" target="_blank">Working example</a>


<b>Chi-square test for independence</b>

`bioinfokit.analys.chisq(table)`

Parameters | Description
------------ | -------------
`table` | CSV delimited data file. It should be contingency table.

Returns:

summary output and variable mosaic plot (mosaic.png)

<a href="https://reneshbedre.github.io/blog/chisq.html" target="_blank">Working example</a>


<b>File format conversions</b>

`bioinfokit.analys.format`

Function | Parameters | Description
------------|---------- | -------------
`bioinfokit.analys.format.fqtofa(file)` | `FASTQ file` | Convert FASTQ file into FASTA format
`bioinfokit.analys.format.hmmtocsv(file)` | `HMM file` | Convert HMM text output (from HMMER tool) to CSV format
`bioinfokit.analys.format.tabtocsv(file)` | `TAB file` | Convert TAB file to CSV format
`bioinfokit.analys.format.csvtotab(file)` | `CSV file` | Convert CSV file to TAB format


Returns:

Output will be saved in same directory

<a href="https://reneshbedre.github.io/blog/format.html" target="_blank">Working example</a>

<b>One-way ANOVA</b>

`bioinfokit.stat.oanova(table, res, xfac, ph, phalpha)`

Parameters | Description
------------ | -------------
`table` | Pandas dataframe in stacked table format
`res` | Response variable (dependent variable) [string][default: None]
`xfac` | Treatments or groups or factors (independent variable) [string][default: None]
`ph` | perform pairwise comparisons with Tukey HSD test [bool (True or False)] [default: False]
`phalpha` |significance level Tukey HSD test [float (0 to 1)][default: 0.05]


Returns:

ANOVA summary, multiple pairwise comparisons, and assumption tests statistics

<a href="https://reneshbedre.github.io/blog/oanova.html" target="_blank">Working example</a>


<b>Manhatten plot</b>

`bioinfokit.visuz.marker.mhat(df, chr, pv, color, dim, r, ar, gwas_sign_line, gwasp, dotsize, markeridcol, markernames, gfont, valpha)`

Parameters | Description
------------ | -------------
`df` |Pandas dataframe object with atleast SNP, chromosome, and P-values columns
`chr` | Name of a column having chromosome numbers [string][default:None]
`pv` | Name of a column having P-values. Must be numeric column [string][default:None]
`color` | List the name of the colors to be plotted. It can accept two alternate colors or the number colors equal to chromosome number. If nothing (None) provided, it will randomly assign the color to each chromosome [list][default:None]
`dim` | Figure size [tuple of two floats (width, height) in inches][default: (6, 4)]
`r` | Figure resolution in dpi [int][default: 300]
`ar` | Rotation of X-axis labels [float][default: 90]
`gwas_sign_line` |Plot statistical significant threshold line defined by option `gwasp` [bool (True or False)][default: False]
`gwasp` |  Statistical significant threshold to identify significant SNPs [float][default: 5E-08]
`dotsize`| The size of the dots in the plot [float][default: 8]
`markeridcol` | Name of a column having SNPs. This is necessary for plotting SNP names on the plot [string][default: None]
`markernames` | The list of the SNPs to display on the plot. These SNP should be present in SNP column. Additionally, it also accepts the dict of SNPs and its associated gene name. If this option set to True, it will label all SNPs with P-value significant score defined by `gwasp` [string, list, dict][default: True]
`gfont` | Font size for SNP names to display on the plot [float][default: 8]
`valpha` | Transparency of points on plot [float (between 0 and 1)][default: 1.0]


Returns:

Manhatten plot image in same directory (manhatten.png)

<a href="https://reneshbedre.github.io/blog/manhat.html" target="_blank">Working example</a>


<b>Extract the sequences from the FASTA file</b>

`bioinfokit.analys.extract_seq(file, id)`

Parameters | Description
------------ | -------------
`file` | input FASTA file from which sequneces to be extracted
`id` | sequence ID file

Returns:
Extracted sequences in FASTA format file in same directory (out.fasta)

<!--
<a href="https://reneshbedre.github.io/blog/extrseq.html" target="_blank">Example</a>
-->


<b>Bar-dot plot</b>

`bioinfokit.visuz.stat.bardot(df, colorbar, colordot, bw, dim, r, ar, hbsize, errorbar, dotsize, markerdot, valphabar, valphadot)`

Parameters | Description
------------ | -------------
`df` |Pandas dataframe object
`colorbar` | Color of bar graph [string or list][default:"#bbcfff"]
`colordot` | Color of dots on bar [string or list][default:"#ee8972"]
`bw` |Width of bar [float][default: 0.4]
`dim` | Figure size [tuple of two floats (width, height) in inches][default: (6, 4)]
`r` | Figure resolution in dpi [int][default: 300]
`ar` | Rotation of X-axis labels [float][default: 0]
`hbsize` | Horizontal bar size for standard error bars [float][default: 4]
`errorbar` |  Draw standard error bars [bool (True or False)][default: True]
`dotsize`| The size of the dots in the plot [float][default: 6]
`markerdot` | Shape of the dot marker. See more options at  https://matplotlib.org/3.1.1/api/markers_api.html [string][default: "o"]
`valphabar` | Transparency of bars on plot [float (between 0 and 1)][default: 1]
`valphadot` | Transparency of dots on plot [float (between 0 and 1)][default: 1]


Returns:

Bra-dot plot image in same directory (bardot.png)

