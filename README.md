
## ChIPexoQual

The idea of this codes is to build a package with a quality control
pipeline for ChIP-exo data sets.

The pipeline considers:

1 - Given a set of reads (in bam file format), is going to build a
collection of regions that separates the genome

2 - For each region is going to calculate a summary statistic designed
to asses a particular bias in the generation of ChIP-exo data set

3 - With the summary statistics is going to build diagnostic plots
that asses the quality of the data set

The package contains other tools that aren't specific to the pipeline,
but useful to asses whether the pipeline is effective as an quality
assessment tool like:

1 - Label the regions in the partition when they overlap any of the
peaks in a previously defined peak list from a ChIP-seq experiment
(with the same conditions)

2 - Calculate some of the current ENCODE quality metrics for a
ChIP-exp dataset


