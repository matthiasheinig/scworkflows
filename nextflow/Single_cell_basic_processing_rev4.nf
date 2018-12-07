params.outdir = 'results'  


if (!params.maxJobs){params.maxJobs = ""} 
if (!params.ncpus){params.ncpus = ""} 
if (!params.samplesheet){params.samplesheet = ""} 

Channel.value(params.maxJobs).set{g_5_maxjobs_g_27}
Channel.value(params.ncpus).into{g_10_ncpus_g_27;g_10_ncpus_g_29}
g_15_xlsfile_g_22 = file(params.samplesheet) 


process split_sample_table_cr_pipe {

publishDir params.outdir, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.txt$/) "samplefiles/$filename"
}

input:
 file sampletable from g_15_xlsfile_g_22

output:
 file '*.txt' into g_22_txtfile_g_20

"""
#!/usr/bin/env Rscript

# Read the table
sampletable="${sampletable}"
table = read.table(sampletable, stringsAsFactors = F, header = T, sep="\t")
for(i in 1:nrow(table)){
  
# Parse params
samplename = table[i,"Sample"]
ncells =  table[i,"ncells"]
fastqpath = table[i,"fastqpath"]
refindex = table[i,"refindex"]
refgtf = table[i,"refgtf"]
chemistry = table[i,"chemistry"]

# fn
fn=paste0(samplename,".txt")

# write to file
fileConn<-file(fn)
writeLines(paste(samplename,ncells,fastqpath,refindex,refgtf,chemistry,sep="\t"), fileConn)
close(fileConn)  
}

"""
}


process create_paramset {

input:
 file input from g_22_txtfile_g_20.flatten()

output:
 set val("${input.baseName}"),file(input) into g_20_paramset_g_27

"""

"""
}


process cellranger_count {

publishDir params.outdir, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /${samplename}$/) "cellranger/$filename"
}

input:
 set samplename,paramfile from g_20_paramset_g_27
 val maxjobs from g_5_maxjobs_g_27
 val samtoolsncpus from g_10_ncpus_g_27

output:
 file "${samplename}" into g_27_folder_g_28, g_27_folder_g_29

"""
# Read paramfile into a var
paramline=\$(<${paramfile})

# Split line
IFS=\$'\t'
params=(\${paramline})

# Create param vars
samplename=\${params[0]}
ncells=\${params[1]}
fastqpath=\${params[2]}
refindex=\${params[3]}
chemistry=\${params[5]}

# Get cellranger version
cellranger sitecheck > sitecheck.txt
crversion=`head -n2 sitecheck.txt | grep -oP '\\(\\K[^)]+'`

cr_cmd="cellranger count --id=cellranger --fastqs=\${fastqpath} \
--sample=${samplename} \
--transcriptome=\${refindex} \
--expect-cells=\${ncells} \
--jobmode=sge --maxjobs=${maxjobs} \
--chemistry=\${chemistry}"

# Echo the cmd to a file
echo \${cr_cmd} > cr_cmd.txt

# Execute
eval \${cr_cmd}

# Create output directories
mkdir -p ${samplename}/count_matrices
mkdir -p ${samplename}/bamfile
mkdir -p ${samplename}/statsfiles

# Copy count matrices folders and h5 files from cellranger directory
cp -R cellranger/outs/filtered_gene_bc_matrices ${samplename}/count_matrices
cp -R cellranger/outs/raw_gene_bc_matrices ${samplename}/count_matrices
cp cellranger/outs/filtered_gene_bc_matrices_h5.h5 ${samplename}/count_matrices
cp cellranger/outs/raw_gene_bc_matrices_h5.h5 ${samplename}/count_matrices

# Copy bamfile
cp cellranger/outs/possorted_genome_bam.bam ${samplename}/bamfile
cp cellranger/outs/possorted_genome_bam.bam.bai ${samplename}/bamfile

# Copy statsfiles
cp cellranger/outs/web_summary.html ${samplename}/statsfiles
cp cellranger/outs/metrics_summary.csv ${samplename}/statsfiles

# Copy analysis folder
cp -R cellranger/outs/analysis ${samplename}

# Copy molecule_info.h5 file
cp cellranger/outs/molecule_info.h5 ${samplename}/analysis

# Try to copy cloupe.cloupe file to analysis folder / is not created if mult genomes are used / do not err if file does not exist
cp cellranger/outs/cloupe.cloupe ${samplename}/analysis 2>/dev/null || :

# Echo the version to a file
echo \${crversion} > ${samplename}/cr_version.txt
# Cp sitecheck file
cp sitecheck.txt ${samplename}/cr_sitecheck.txt
# Copy the results directory
cp cr_cmd.txt ${samplename}/cr_cmd.txt

# Get samtools version
samtools --help | head -3 > ${samplename}/samtools_version.txt

# Sort bamfile by CB tag / for velocyto
samt_cmd="samtools sort --threads ${samtoolsncpus} -t CB -O BAM -o ${samplename}/bamfile/cellsorted_possorted_genome_bam.bam ${samplename}/bamfile/possorted_genome_bam.bam"

# Execute
eval \${samt_cmd}

# Echo cmd to file
echo \${samt_cmd} > ${samplename}/samtools_cmd.txt

# Copy the param file to the out directory
cp ${paramfile} ${samplename}/params.txt

# rm cellranger folder
#rm -R cellranger

"""
}


process Velocyto_folder {

publishDir params.outdir, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /${name}$/) "velocyto/$filename"
}

input:
 file name from g_27_folder_g_29
 val ncpus from g_10_ncpus_g_29

output:
 file "${name}" into g_29_folder

"""
# Read the paramsfile and split to get gft file

export LC_ALL=C.UTF-8
export LANG=C.UTF-8


# Read paramfile into a var
paramline=\$(<${name}/params.txt)

# Get org IFS
OLDIFS=\$IFS

# Split line
IFS=\$'\t'
params=(\${paramline})

# Reset IFS
IFS=\$OLDIFS

# gtf is 5. param (note 0 based indexing)
gtffile=\${params[4]}

# Get version
velocyto --version > version.txt
vcversion=`egrep -o "([0-9]{1,}\\.)+[0-9]{1,}" version.txt`

# run velocyto for each genome subdirectory
# Generate an array of subdirs in folder
dirarray=( `ls ${name}/count_matrices/filtered_gene_bc_matrices/` )

# Iterate over the folders
for dname in "\${dirarray[@]}"
do
# Get barcode file path
bcpath="${name}/count_matrices/filtered_gene_bc_matrices/\${dname}/barcodes.tsv"

# Run velocyto
velocmd="velocyto run --samtools-threads ${ncpus} -b \${bcpath} -o velocyto/\${dname} ${name}/bamfile/possorted_genome_bam.bam \${gtffile}"

# create results dir
mkdir -p "velocyto/\${dname}"

# Echo cmd to file
echo \${velocmd} > "velocyto/\${dname}/vc.cmd.txt"

# Execute
eval \${velocmd}
done

# mv the original folder (which is a softlink)
mv ${name} oldname

# mv velocyto folder to name
mv velocyto ${name}

# cp velocyto version file
cp version.txt ${name}/vc_version.txt


"""
}


process Create_AnnData {

publishDir params.outdir, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /${name}$/) "scanpy_AnnData/$filename"
}

input:
 file name from g_27_folder_g_28

output:
 file "${name}" into g_28_folder

"""
#!/usr/bin/env python
# Adapted from https://nbviewer.jupyter.org/github/theislab/scanpy_usage/blob/master/170505_seurat/seurat.ipynb

import numpy as np
import pandas as pd
import scanpy.api as sc
import glob
import os
import sys

# rename directory
os.rename('${name}',"cellranger")


# Create annData for filtered results
# Glob genome folders
pathg =glob.glob('./cellranger/count_matrices/filtered_gene_bc_matrices/*')

# Iterate over the genome folder
for x in pathg:
# Get the basename
    genomeFname=os.path.basename(x)

    adata = sc.read('./cellranger/count_matrices/filtered_gene_bc_matrices/' + genomeFname + '/matrix.mtx', cache=True).T  # transpose the data
    adata.var_names = pd.read_csv('./cellranger/count_matrices/filtered_gene_bc_matrices/' + genomeFname + '/genes.tsv', header=None, sep='\t')[1]
    adata.obs_names = pd.read_csv('./cellranger/count_matrices/filtered_gene_bc_matrices/' + genomeFname + '/barcodes.tsv', header=None)[0]
    adata.var_names_make_unique()

    # Write to file
    results_file = './${name}/' + genomeFname + '/filtered_gene_bc_matrices.h5ad'
    adata.write(results_file)

    # Create annData for full results
    # pathg =glob.glob('./cellranger/count_matrices/raw_gene_bc_matrices/*/')
    # path=str(pathg[0])

    adata = sc.read('./cellranger/count_matrices/raw_gene_bc_matrices/' + genomeFname + '/matrix.mtx', cache=True).T  # transpose the data
    adata.var_names = pd.read_csv('./cellranger/count_matrices/raw_gene_bc_matrices/' + genomeFname + '/genes.tsv', header=None, sep='\t')[1]
    adata.obs_names = pd.read_csv('./cellranger/count_matrices/raw_gene_bc_matrices/' + genomeFname + '/barcodes.tsv', header=None)[0]
    adata.var_names_make_unique()

    # Write to file
    results_file = './${name}/' + genomeFname + '/raw_gene_bc_matrices.h5ad'
    adata.write(results_file)

# Write version of sc to file
import sys, io
stdout = sys.stdout
sys.stdout = io.StringIO()
sc.logging.print_versions()
# get output and restore sys.stdout
output = sys.stdout.getvalue()
sys.stdout = stdout
# write output to file
file = open('./${name}/scanpy_version.txt',"w") 
file.write(output)
"""
}


workflow.onComplete {
println "##Pipeline execution summary##"
println "---------------------------"
println "##Completed at: $workflow.complete"
println "##Duration: ${workflow.duration}"
println "##Success: ${workflow.success ? 'OK' : 'failed' }"
println "##Exit status: ${workflow.exitStatus}"
}
