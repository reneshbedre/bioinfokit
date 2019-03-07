Bioinformatics data analysis and visualization toolkit


**<span style="color:#33a8ff">How to install:</span>**

For github,

```
git clone https://github.com/reneshbedre/bioinfokit.git
cd bioinfokit
python3 setup.py install
```

From Python package index,

`pip install -i https://test.pypi.org/simple/ bioinfokit `

**<span style="color:#33a8ff">Functions:</span>**

<b>Volcano plot</b>

`bioinfokit.visuz.volcano(table, lfc, pv, lfc_thr, pv_thr)`

Parameters | Description
------------ | -------------
`table` |Comma separated (csv) gene expression table having atleast gene IDs, log fold change, P-values or adjusted P-values
`lfc` | Name of a column having log fold change values [default:logFC]
`pv` | Name of a column having P-values or adjusted P-values [default:p_values]
`lfc_thr` | Log fold change cutoff for up and downregulated genes [default:1]
`pv_thr` | P-values or adjusted P-values cutoff for up and downregulated genes [default:0.05]

Returns:

Volcano plot image in same directory (volcano.png)

Example:
```
Python 3.6.4 (default, Oct 10 2018, 14:08:09)
[GCC Intel(R) C++ gcc 6.4 mode] on linux
Type "help", "copyright", "credits" or "license" for more information.
>>> from bioinfokit import visuz

```


<b>MA plot</b>

`bioinfokit.visuz.ma(table, lfc, ct_count, st_count, pv_thr)`

Parameters | Description
------------ | -------------
`table` |Comma separated (csv) gene expression table having atleast gene IDs, log fold change, P-values or adjusted P-values
`lfc` | Name of a column having log fold change values [default:logFC]
`ct_count` | Name of a column having count values for control sample [default:value1]
`st_count` | Name of a column having count values for treatment sample [default:value2]
`lfc_thr` | Log fold change cutoff for up and downregulated genes [default:1]

Returns:

MA plot image in same directory (ma.png)


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




<b>Correlation matrix plot</b>

`bioinfokit.visuz.involcano(table, lfc, pv, lfc_thr, pv_thr)`

Parameters | Description
------------ | -------------
`table` |Comma separated (csv) table with numerical variables (columns) to find correlation. Ideally, you should have three or more variables. Table should not have identifier column.
`corm` | Correlation method [pearson,kendall,spearman] [default:pearson]


Returns:

Correlation matrix plot image in same directory (corr_mat.png)