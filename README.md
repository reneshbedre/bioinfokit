
[![DOI](https://zenodo.org/badge/174428856.svg)](https://zenodo.org/badge/latestdoi/174428856)


The bioinfokit toolkit aimed to provide various easy-to-use functionalities to analyze,  
visualize, and interpret the biological data generated from genome-scale omics experiments.

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
python setup.py install
```

<!--
From Python package index,

`pip install -i https://test.pypi.org/simple/ bioinfokit `

**<span style="color:#33a8ff">Functions:</span>**
-->


<b>Volcano plot</b>

`bioinfokit.visuz.gene_exp.volcano(table, lfc, pv, lfc_thr, pv_thr, color, valpha, geneid, genenames, gfont, gstyle, sign_line,
    dotsize, markerdot, r, dim, show, figtype, axxlabel, axylabel, axlabelfontsize, axlabelfontname, axtickfontsize, axtickfontname,
    xlm, ylm)`

Parameters | Description
------------ | -------------
`table` |Pandas dataframe table having atleast gene IDs, log fold change, P-values or adjusted P-values columns
`lfc` | Name of a column having log fold change values [string][default:logFC]
`pv` | Name of a column having P-values or adjusted P-values [string][default:p_values]
`lfc_thr` | Log fold change cutoff for up and downregulated genes [float][default:1.0]
`pv_thr` | P-values or adjusted P-values cutoff for up and downregulated genes [float][default:0.05]
`color` | Tuple of two colors [tuple][default: ("green", "red")]
`valpha` | Transparency of points on volcano plot [float (between 0 and 1)][default: 1.0]
`geneid` | Name of a column having gene Ids. This is necessary for plotting gene label on the points [string][default: None]
`genenames` | Tuple of gene Ids to label the points. The gene Ids must be present in the geneid column. If this option set to "deg" it will label all genes defined by lfc_thr and pv_thr [string, tuple, dict][default: None]
`gfont` | Font size for genenames [float][default: 10.0]. gfont not compatible with gstyle=2.
`gstyle` | Style of the text for genenames. 1 for default text and 2 for box text [int][default: 1]
`sign_line` | Show grid lines on plot with defined log fold change (`lfc_thr`) and P-value (`pv_thr`) threshold value [True or False][default:False]
`dotsize`| The size of the dots in the plot [float][default: 8]
`markerdot` | Shape of the dot marker. See more options at  https://matplotlib.org/3.1.1/api/markers_api.html [string][default: "o"]
`dim` | Figure size [tuple of two floats (width, height) in inches][default: (5, 5)]
`r` | Figure resolution in dpi [int][default: 300]. Not compatible with `show`= True
`figtype` | Format of figure to save. Supported format are eps, pdf, pgf, png, ps, raw, rgba, svg, svgz [string][default:'png']
`show`  | Show the figure on console instead of saving in current folder [True or False][default:False]
`axxlabel` | Label for X-axis. If you provide this option, default label will be replaced [string][default: None]
`axylabel` | Label for Y-axis. If you provide this option, default label will be replaced [string][default: None]
`axlabelfontsize` | Font size for axis labels [float][default: 9]
`axlabelfontname` | Font name for axis labels [string][default: 'Arial']
`axtickfontsize` | Font size for axis ticks [float][default: 9]
`axtickfontname` | Font name for axis ticks [string][default: 'Arial']
`xlm` | Range of ticks to plot on X-axis [float (left, right, interval)][default: None]
`ylm` | Range of ticks to plot on Y-axis [float (bottom, top, interval)][default: None]

Returns:

Volcano plot image in same directory (volcano.png)
<a href="https://reneshbedre.github.io/blog/volcano.html" target="_blank">Working example</a>

<b>MA plot</b>

`bioinfokit.visuz.gene_exp.ma(table, lfc, ct_count, st_count, lfc_thr, color, dim, dotsize, show, r, valpha, figtype, axxlabel,
    axylabel, axlabelfontsize, axtickfontsize, axtickfontname, xlm, ylm)`

Parameters | Description
------------ | -------------
`table` |Pandas dataframe  table having atleast gene IDs, log fold change, and counts (control and treatment) columns
`lfc` | Name of a column having log fold change values [default:logFC]
`ct_count` | Name of a column having count values for control sample [default:value1]
`st_count` | Name of a column having count values for treatment sample [default:value2]
`lfc_thr` | Log fold change cutoff for up and downregulated genes [default:1]
`color` | Tuple of two colors [tuple][default: ("green", "red")]
`dotsize`| The size of the dots in the plot [float][default: 8]
`markerdot` | Shape of the dot marker. See more options at  https://matplotlib.org/3.1.1/api/markers_api.html [string][default: "o"]
`valpha` | Transparency of points on plot [float (between 0 and 1)][default: 1.0]
`dim` | Figure size [tuple of two floats (width, height) in inches][default: (5, 5)]
`r` | Figure resolution in dpi [int][default: 300]. Not compatible with `show`= True
`figtype` | Format of figure to save. Supported format are eps, pdf, pgf, png, ps, raw, rgba, svg, svgz [string][default:'png']
`show`  | Show the figure on console instead of saving in current folder [True or False][default:False]
`axxlabel` | Label for X-axis. If you provide this option, default label will be replaced [string][default: None]
`axylabel` | Label for Y-axis. If you provide this option, default label will be replaced [string][default: None]
`axlabelfontsize` | Font size for axis labels [float][default: 9]
`axtickfontsize` | Font size for axis ticks [float][default: 9]
`axtickfontname` | Font name for axis ticks [string][default: 'Arial']
`xlm` | Range of ticks to plot on X-axis [float (left, right, interval)][default: None]
`ylm` | Range of ticks to plot on Y-axis [float (bottom, top, interval)][default: None]

Returns:

MA plot image in same directory (ma.png)
<a href="https://reneshbedre.github.io/blog/ma.html" target="_blank">Working example</a>

<b>Inverted Volcano plot</b>

`bioinfokit.visuz.gene_exp.involcano(table, lfc, pv, lfc_thr, pv_thr, color, valpha, geneid, genenames, gfont, gstyle,
    dotsize, markerdot, r, dim, show, r, dim, show, figtype, axxlabel, axylabel, axlabelfontsize, axtickfontsize, axtickfontname)`

Parameters | Description
------------ | -------------
`table` |Pandas dataframe table having atleast gene IDs, log fold change, P-values or adjusted P-values
`lfc` | Name of a column having log fold change values [default:logFC]
`pv` | Name of a column having P-values or adjusted P-values [default:p_values]
`lfc_thr` | Log fold change cutoff for up and downregulated genes [default:1]
`pv_thr` | P-values or adjusted P-values cutoff for up and downregulated genes [default:0.05]
`color` | Tuple of two colors [tuple][default: ("green", "red")]
`valpha` | Transparency of points on volcano plot [float (between 0 and 1)][default: 1.0]
`geneid` | Name of a column having gene Ids. This is necessary for plotting gene label on the points [string][default: None]
`genenames` | Tuple of gene Ids to label the points. The gene Ids must be present in the geneid column. If this option set to "deg" it will label all genes defined by lfc_thr and pv_thr [string, tuple, dict][default: None]
`gfont` | Font size for genenames [float][default: 10.0]
`gstyle` | Style of the text for genenames. 1 for default text and 2 for box text [int][default: 1]
`dotsize`| The size of the dots in the plot [float][default: 8]
`markerdot` | Shape of the dot marker. See more options at  https://matplotlib.org/3.1.1/api/markers_api.html [string][default: "o"]
`dim` | Figure size [tuple of two floats (width, height) in inches][default: (5, 5)]
`r` | Figure resolution in dpi [int][default: 300]. Not compatible with `show`= True
`figtype` | Format of figure to save. Supported format are eps, pdf, pgf, png, ps, raw, rgba, svg, svgz [string][default:'png']
`show`  | Show the figure on console instead of saving in current folder [True or False][default:False]
`axxlabel` | Label for X-axis. If you provide this option, default label will be replaced [string][default: None]
`axylabel` | Label for Y-axis. If you provide this option, default label will be replaced [string][default: None]
`axlabelfontsize` | Font size for axis labels [float][default: 9]
`axtickfontsize` | Font size for axis ticks [float][default: 9]
`axtickfontname` | Font name for axis ticks [string][default: 'Arial']


Returns:

Inverted volcano plot image in same directory (involcano.png)

<a href="https://reneshbedre.github.io/blog/volcano.html" target="_blank">Working example</a>


<b>Correlation matrix plot</b>

`bioinfokit.visuz.stat.corr_mat(table, corm, cmap, r, dim, show, figtype, axtickfontsize, axtickfontname)`

Parameters | Description
------------ | -------------
`table` | Dataframe object with numerical variables (columns) to find correlation. Ideally, you should have three or more variables. Dataframe should not have identifier column.
`corm` | Correlation method [pearson,kendall,spearman] [default:pearson]
`cmap` | Color Palette for heatmap [string][default: 'seismic']. More colormaps are available at  
         https://matplotlib.org/3.1.0/tutorials/colors/colormaps.html
`r` | Figure resolution in dpi [int][default: 300]. Not compatible with `show`= True
`dim` | Figure size [tuple of two floats (width, height) in inches][default: (6, 5)]        
`show`  | Show the figure on console instead of saving in current folder [True or False][default:False]
`figtype` | Format of figure to save. Supported format are eps, pdf, pgf, png, ps, raw, rgba, svg, svgz [string][default:'png']
`axtickfontsize` | Font size for axis ticks [float][default: 7]
`axtickfontname` | Font name for axis ticks [string][default: 'Arial']

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

`bioinfokit.visuz.gene_exp.hmap(table, cmap='seismic', scale=True, dim=(6, 8), clus=True, zscore=None, xlabel=True, ylabel=True, 
    tickfont=(12, 12), show, r, figtype)`

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
`show`  | Show the figure on console instead of saving in current folder [True or False][default:False]
`r` | Figure resolution in dpi [int][default: 300]. Not compatible with `show`= True
`figtype` | Format of figure to save. Supported format are eps, pdf, pgf, png, ps, raw, rgba, svg, svgz [string][default:'png']

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

<b>Two sample and Welch's t-test </b>

`bioinfokit.analys.stat.ttsam(table, xfac, res, evar, alpha)`

Parameters | Description
------------ | -------------
`table` | Pandas dataframe. It should be stacked table with independent (xfac) and dependent (res) variable columns.
`xfac` | Independent group column name with two levels [string][default: None]
`res` | Response variable column name [string][default: None]
`evar` | t-test with equal variance [bool (True or False)][default: True]
`alpha` | Confidence level [float][default: 0.05]

Returns:

summary output and group boxplot (ttsam_boxplot.png)

<a href="https://reneshbedre.github.io/blog/ttest.html" target="_blank">Working example</a>


<b>Chi-square test for independence</b>

`bioinfokit.analys.stat.chisq(table)`

Parameters | Description
------------ | -------------
`table` | Pandas dataframe. It should be contingency table.

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

`bioinfokit.visuz.marker.mhat(df, chr, pv, color, dim, r, ar, gwas_sign_line, gwasp, dotsize, markeridcol, markernames, 
    gfont, valpha, show, figtype, axxlabel, axylabel, axlabelfontsize, ylm)`

Parameters | Description
------------ | -------------
`df` |Pandas dataframe object with atleast SNP, chromosome, and P-values columns
`chr` | Name of a column having chromosome numbers [string][default:None]
`pv` | Name of a column having P-values. Must be numeric column [string][default:None]
`color` | List the name of the colors to be plotted. It can accept two alternate colors or the number colors equal to chromosome number. If nothing (None) provided, it will randomly assign the color to each chromosome [list][default:None]
`gwas_sign_line` |Plot statistical significant threshold line defined by option `gwasp` [bool (True or False)][default: False]
`gwasp` |  Statistical significant threshold to identify significant SNPs [float][default: 5E-08]
`dotsize`| The size of the dots in the plot [float][default: 8]
`markeridcol` | Name of a column having SNPs. This is necessary for plotting SNP names on the plot [string][default: None]
`markernames` | The list of the SNPs to display on the plot. These SNP should be present in SNP column. Additionally, it also accepts the dict of SNPs and its associated gene name. If this option set to True, it will label all SNPs with P-value significant score defined by `gwasp` [string, list, dict][default: True]
`gfont` | Font size for SNP names to display on the plot [float][default: 8]
`valpha` | Transparency of points on plot [float (between 0 and 1)][default: 1.0]
`dim` | Figure size [tuple of two floats (width, height) in inches][default: (6, 4)]
`r` | Figure resolution in dpi [int][default: 300]
`ar` | Rotation of X-axis labels [float][default: 90]
`figtype` | Format of figure to save. Supported format are eps, pdf, pgf, png, ps, raw, rgba, svg, svgz [string][default:'png']
`show`  | Show the figure on console instead of saving in current folder [True or False][default:False]
`axxlabel` | Label for X-axis. If you provide this option, default label will be replaced [string][default: None]
`axylabel` | Label for Y-axis. If you provide this option, default label will be replaced [string][default: None]
`axlabelfontsize` | Font size for axis labels [float][default: 9]
`ylm` | Range of ticks to plot on Y-axis [float tuple (bottom, top, interval)][default: None]


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

<b>Bar-dot plot</b>

`bioinfokit.visuz.stat.bardot(df, colorbar, colordot, bw, dim, r, ar, hbsize, errorbar, dotsize, markerdot, valphabar, 
    valphadot, show, figtype, axxlabel, axylabel, axlabelfontsize, axlabelfontname, ylm, axtickfontsize, axtickfontname,
    yerrlw, yerrcw)`

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
`figtype` | Format of figure to save. Supported format are eps, pdf, pgf, png, ps, raw, rgba, svg, svgz [string][default:'png']
`show`  | Show the figure on console instead of saving in current folder [True or False][default:False]
`axxlabel` | Label for X-axis. If you provide this option, default label will be replaced [string][default: None]
`axylabel` | Label for Y-axis. If you provide this option, default label will be replaced [string][default: None]
`axlabelfontsize` | Font size for axis labels [float][default: 9]
`axlabelfontname` | Font name for axis labels [string][default: 'Arial']
`ylm` | Range of ticks to plot on Y-axis [float tuple (bottom, top, interval)][default: None]
`axtickfontsize` | Font size for axis ticks [float][default: 9]
`axtickfontname` | Font name for axis ticks [string][default: 'Arial']
`yerrlw` | Error bar line width [float][default: None]
`yerrcw` | Error bar cap width [float][default: None]


Returns:

Bar-dot plot image in same directory (bardot.png)

<a href="https://reneshbedre.github.io/blog/bardot.html" target="_blank">Working Example</a>


<b>FASTQ quality format detection</b>

`bioinfokit.analys.format.fq_qual_var(file)`

Parameters | Description
------------ | -------------
`file` |FASTQ file to detect quality format [deafult: None]

Returns:

Quality format encoding name for FASTQ file (Supports only Sanger, Illumina 1.8+ and Illumina  1.3/1.4)

<a href="https://reneshbedre.github.io/blog/fqqualfmt.html" target="_blank">Working Example</a>


<b>Linear regression analysis</b>

`bioinfokit.visuz.stat.lin_reg(df, x, y)`

Parameters | Description
------------ | -------------
`df` |Pandas dataframe object
`x` | Name of column having independent X variables [list][default:None]
`y` | Name of column having dependent Y variables [list][default:None]

Returns:

Regression analysis summary

<a href="https://reneshbedre.github.io/blog/linearreg.html" target="_blank">Working Example</a>


<b>Regression plot</b>

`bioinfokit.visuz.stat.regplot(df, x, y, yhat, dim, colordot, colorline, r, ar, dotsize, markerdot, linewidth, 
    valphaline, valphadot, show, figtype, axxlabel, axylabel, axlabelfontsize, axlabelfontname, xlm, ylm, axtickfontsize,
    axtickfontname)`

Parameters | Description
------------ | -------------
`df` | Pandas dataframe object
`x` | Name of column having independent X variables [string][default:None]
`y` | Name of column having dependent Y variables [string][default:None]
`yhat` |Name of column having predicted response of Y variable (y_hat) from regression [string][default:None]
`dim` | Figure size [tuple of two floats (width, height) in inches][default: (6, 4)]
`r` | Figure resolution in dpi [int][default: 300]
`ar` | Rotation of X-axis labels [float][default: 0]
`dotsize`| The size of the dots in the plot [float][default: 6]
`markerdot` | Shape of the dot marker. See more options at  https://matplotlib.org/3.1.1/api/markers_api.html [string][default: "o"]
`valphaline` | Transparency of regression line on plot [float (between 0 and 1)][default: 1]
`valphadot` | Transparency of dots on plot [float (between 0 and 1)][default: 1]
`linewidth` | Width of regression line [float][default: 1]
`figtype` | Format of figure to save. Supported format are eps, pdf, pgf, png, ps, raw, rgba, svg, svgz [string][default:'png']
`show`  | Show the figure on console instead of saving in current folder [True or False][default:False]
`axxlabel` | Label for X-axis. If you provide this option, default label will be replaced [string][default: None]
`axylabel` | Label for Y-axis. If you provide this option, default label will be replaced [string][default: None]
`axlabelfontsize` | Font size for axis labels [float][default: 9]
`axlabelfontname` | Font name for axis labels [string][default: 'Arial']
`xlm` | Range of ticks to plot on X-axis [float tuple (bottom, top, interval)][default: None]
`ylm` | Range of ticks to plot on Y-axis [float tuple (bottom, top, interval)][default: None]
`axtickfontsize` | Font size for axis ticks [float][default: 9]
`axtickfontname` | Font name for axis ticks [string][default: 'Arial']


Returns:

Regression plot image in same directory (reg_plot.png)

<a href="https://reneshbedre.github.io/blog/linearreg.html" target="_blank">Working Example</a>




References:
- Travis E. Oliphant. A guide to NumPy, USA: Trelgol Publishing, (2006).
- John D. Hunter. Matplotlib: A 2D Graphics Environment, Computing in Science & Engineering, 9, 90-95 (2007), DOI:10.1109/MCSE.2007.55 (publisher link)
- Fernando Pérez and Brian E. Granger. IPython: A System for Interactive Scientific Computing, Computing in Science & Engineering, 9, 21-29 (2007), DOI:10.1109/MCSE.2007.53 (publisher link)
- Michael Waskom, Olga Botvinnik, Joel Ostblom, Saulius Lukauskas, Paul Hobson, MaozGelbart, … Constantine Evans. (2020, January 24). mwaskom/seaborn: v0.10.0 (January 2020) (Version v0.10.0). Zenodo. http://doi.org/10.5281/zenodo.3629446
- Fabian Pedregosa, Gaël Varoquaux, Alexandre Gramfort, Vincent Michel, Bertrand Thirion, Olivier Grisel, Mathieu Blondel, Peter Prettenhofer, Ron Weiss, Vincent Dubourg, Jake Vanderplas, Alexandre Passos, David Cournapeau, Matthieu Brucher, Matthieu Perrot, Édouard Duchesnay. Scikit-learn: Machine Learning in Python, Journal of Machine Learning Research, 12, 2825-2830 (2011)
- Wes McKinney. Data Structures for Statistical Computing in Python, Proceedings of the 9th Python in Science Conference, 51-56 (2010)

bioinfokit cited by:
- Greaney AM, Adams TS, Raredon MS, Gubbins E, Schupp JC, Engler AJ, Ghaedi M, Yuan Y, Kaminski N, Niklason LE. Platform 
  Effects on Regeneration by Pulmonary Basal Cells as Evaluated by Single-Cell RNA Sequencing. Cell Reports. 2020 Mar 
  24;30(12):4250-65.