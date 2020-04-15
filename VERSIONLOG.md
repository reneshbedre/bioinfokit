v0.7 has the following updates and changes
- `split_fastq` function added for splitting individual (left and right) paired-end fastq files
  from single interleaved paired-end file
- GFF3 to GTF file conversion utility added under class `gff`
- two-sample t-test updated for CI
- module termcolor removed

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
