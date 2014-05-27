# What can you do with R
The tasks that are descrived above can be accomplished by R. R is not only a powerful statistical programming language but also go-to data analysis tool for many computational genomics experts. High-dimensional genomics datasets are usually suitable to be analyzed with core R packages and functions. On top of that, Bioconductor and CRAN have an array of specialized tools for doing genomics specific analysis.
Here is a list of computational genomics tasks that can be completed using R.

### Data munging (pre-processing)

Often times, the data will not come in ready to analyze format. You may need to convert it to other formats by transforming data points (such as log transforming, normalizing etc), or remove columns/rows that , or remove data points with empty values and, and subset the data set with some arbitrary condition. Most of these tasks can be achieved using R. In addition, with the help of packages R can con- nect to databases in various formats such as mySQL, mongoDB, etc., and query and get the data to R environment using database specific tools. Unfortunately, not all data muging and processing tasks can
be accomplished only by R. At times, you may need to use domain specific software or software dealing better with specific type of data sets. For example, R is not great at dealing with character strings, if you are trying to filter a large dataset based on some regular expres- sion you may be better of with perl or awk.

### General data anaylsis and exploration
Most genomics data sets are suitable for application of general data analysis tools. In some cases, you may need to preprocess the data to get it to a state that is suitable for application such tools.
* unsupervised data analysis: clustering (k-means, hierarchical), matrix factorization (PCA, ICA etc)
* supervised data analysis: generalized linear models, support vec- tor machines, randomForests


### Visualization
Visualization is an important part of all data analysis techniques including computational genomics. Again, you can use core visu- alization technniques in R and also genomics specific ones with the help of specific packages.
* Basic plots: Histograms, scatter plots, bar plots, box plots
* ideograms and circus plots for genomics
* heatmaps
* meta-profiles of genomic features, read enrichment over all pro- moters
* genomic track visualization for given locus


### Dealing with genomic intervals
Most of the genomics data come in a tabular format that contains the location in the genome and some other relevant values, such as scores for those genomic features and/or names. R/Bioconductor has dedicated methods to deal with such data. Here are a couple of example tasks that you can achieve using R.
* Overlapping CpG islands with transcription start sites, and filter- ing based on overlaps
* Aligning reads and making read enrichment profiles
* Overlapping aligned reads with exons and counting aligned reads per gene


### Application of other bioinformatics specific algorithms
In addition to genomic interval centered methods, R/Bioconductor gives you access to multitude of other bioinformatics specific algo- rithms. Here are some of the things you can do.
* Sequence analysis: TF binding motifs, GC content and CpG counts of a given DNA sequence
* Differential expression (or arrays and sequencing based measure- ments)
* Gene set/Pathway analysis: What kind of genes are enriched in my gene set
