# 10X genomics
The 10X genomics pipeline performs basic analysis of data generated using the 10X genomics protocol. It produces gene x cell count matrices using the cell ranger software. In addition an RNA velocity analysis is run via velocyto.

Input files and parameters are summarized in the sample table file, with mandatory fields
* Sample: Name of the sample
* ncells: Estimated number of cells
* fastqpath: Path to directory containing the sequencing reads in fastq format
* refindex: Path to the reference genome index
* refgtf: Filename of a gtf file describing the gene models
* chemistry: Version of the chemistry
