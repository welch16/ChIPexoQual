
ChIPexoQual
===========

ChIPexoQual provides a quality control pipeline of ChIP-exo nexus/data. It allows the quick evaluation of ChIP-exo/nexus data quality by directly operating on aligned read files and without requiring any complex statistical modelling or intensive computation such as identification of potential binding regions/events. This enables its broad application and versatile utility, and provide easy data evaluation before any statistical analysis for identifying binding events. 

The overall pipeline follows the steps:

1. Partitions the reference genome into islands representing overlapping clusters of reads separated by gaps. 
2. Calculates a set of summary statistics.
3. Visualize the data with a collection of diagnostic plots aimed at quantifying ChIP enrichment and strand imbalance.
4. Generates quantitative summaries of these diagnostic plots.

To install the package, it is easiest to use:

```
#install.packages("devtools")
library(devtools)
install_github("welch16/ChIPexoQualExample")
install_github("welch16/ChIPexoQual")
```

References
==========

Welch R, Chung D, Grass J, Landick R, and Keles S. "Data Exploration, Quality Control, and Statistical Analysis of ChIP-exo/nexus Experiments" (in preparation)


