# Dropseq
The dropseq pipeline performs basic analysis of data generated using the dropseq protocol. It produces gene x cell count matrices according to the XX procedure. In addition an RNA velocity analysis is run via velocyto.

Input files and parameters are summarized in the sample table file, with mandatory fields
* SampleName: Name of the sample
* fq1: Filename of the file containing the first reads of a set of read pairs
* fq2: Filename of the file containing second reads of a set of read pairs
* nestcells: Estimated number of cells
* genomePath: Path to the reference genome index
* genomeFasta: Fasta file with the reference genome sequence
* gtf: Filename of a gtf file describing the gene models