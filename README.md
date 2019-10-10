Bioinformatics data analysis and visualization toolkit


**<span style="color:#33a8ff">How to install:</span>**

From GitHub,

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

`bioinfokit.visuz.volcano(table, lfc, pv, lfc_thr, pv_thr)`

Parameters | Description
------------ | -------------
`table` |Comma separated (csv) gene expression table having atleast gene IDs, log fold change, P-values or adjusted P-values columns
`lfc` | Name of a column having log fold change values [default:logFC]
`pv` | Name of a column having P-values or adjusted P-values [default:p_values]
`lfc_thr` | Log fold change cutoff for up and downregulated genes [default:1]
`pv_thr` | P-values or adjusted P-values cutoff for up and downregulated genes [default:0.05]

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

`bioinfokit.visuz.involcano(table, lfc, pv, lfc_thr, pv_thr)`

Parameters | Description
------------ | -------------
`table` |Comma separated (csv) gene expression table having atleast gene IDs, log fold change, P-values or adjusted P-values
`lfc` | Name of a column having log fold change values [default:logFC]
`pv` | Name of a column having P-values or adjusted P-values [default:p_values]
`lfc_thr` | Log fold change cutoff for up and downregulated genes [default:1]
`pv_thr` | P-values or adjusted P-values cutoff for up and downregulated genes [default:0.05]

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

<!--
<a href="https://reneshbedre.github.io/blog/pca_3d.html" target="_blank">Working example</a>

<b>Extract the subsequence from the genome or gene sequences</b>

`bioinfokit.analys.rev_com(sequence)`

Parameters | Description
------------ | -------------
`sequence file` | Genome or gene sequence file in fasta format
`id` | sequence ID
`start` | Start coordinate of the sequence to extract
`end` | End coordinate of the sequence to extract
`strand` | Sequence strand [default: plus]

Returns:
Extracted subsequence

<a href="https://reneshbedre.github.io/blog/pca_3d.html" target="_blank">Working example</a>
-->