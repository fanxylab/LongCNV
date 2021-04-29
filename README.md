# LongCNV

The [**LongCNV**](https://github.com/fanxylab/LongCNV.git) is used for detecting copy number of long reads bam form nanopore and pacbio sequencing.
- We split consecutive bins along the hg38 reference genome with each bin of specified  length after excluding the N-gaps region which downloaded from the [**UCSC**](https://genome.ucsc.edu/cgi-bin/hgTables). 
- Two rounds of locally weighted regression (LOESS) model were used for GC-content correction.The copy number of each bin was defined by the ratio of read counts and predicted values of final model and normalized by log2.
- To infer discrete copy number segments, the circular binary segmentation (CBS) algorithm was adopted to get the mean CN in stable segment region.

### User installation and useage
Refer  [**ecDNAFinder**](https://github.com/fanxylab/ecDNAFinder.git) to get the installation and useage.

