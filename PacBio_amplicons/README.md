# PacBio amplicons
This directory contains the scripts to generate two consensus haplotypes per individual starting from q20 (>99% accuracy) PacBio Circular Consensus Sequences (CCS) for the ACCase and ALS loci. This includes the pacbio amplicon analysis (pbaa) step, and the post-processing of the resulting clusters. The pipeline is exemplified with sample “DE01321_01”.

In addition, this directory contains the custom R-scripts to generate phylogenetic trees of both ACCase and ALS for the 47 populations analyzed. Furthermore at the end of the scripts, a loop to generate the metadata files for the TSR annotation in POPART has been implemented. The necessary multiple alignment input files in nexus format for POPART, and the TBE files inferred by RAXML-NG for the ggtrees are included.
