
## ChIPexoQual

The idea of this codes is to build a package with a quality control
pipeline for ChIP-exo data sets.

The pipeline considers:

1 - Given a set of reads (in bam file format), is going to build a
collection of regions that separates the genome

2 - For each region is going to calculate a collection of summary
statistics to asses biases that are specific to the generation of
ChIP-exo data.

3 - Using those summary statistics it creates diagnostics plots to
help understand the behavior of those biases.



