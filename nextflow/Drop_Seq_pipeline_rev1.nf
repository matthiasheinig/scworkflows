params.outdir = 'results'  


if (!params.samplesheet){params.samplesheet = ""} 
if (!params.ncpu){params.ncpu = ""} 

g_1_xlsfile_g_21 = file(params.samplesheet) 
Channel.value(params.ncpu).set{g_16_ncpus_g_22}


process split_sample_table_dropseq {

publishDir params.outdir, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.txt$/) "samplefiles/$filename"
}

input:
 file sampletable from g_1_xlsfile_g_21

output:
 file '*.txt' into g_21_txtfile_g_7

"""
#!/usr/bin/env Rscript

# Read the table
sampletable="${sampletable}"
table = read.table(sampletable, stringsAsFactors = F, header = T, sep="\t")
for(i in 1:nrow(table)){
  
# Parse params
samplename = table[i,"SampleName"]
fq1 =  table[i,"fq1"]
fq2 =  table[i,"fq2"]
#nestcells = table[i,"nestcells"]
genomepath =  table[i,"genomePath"]
genomefasta =  table[i,"genomeFasta"]
gtf =  table[i,"gtf"]
# fn
fn=paste0(samplename,".txt")

# write to file
fileConn<-file(fn)
writeLines(paste(samplename,fq1,fq2,genomepath,genomefasta,gtf, sep="\t"), fileConn)
close(fileConn)  
}

"""
}


process create_paramset {

input:
 file input from g_21_txtfile_g_7.flatten()

output:
 set val("${input.baseName}"),file(input) into g_7_paramset_g_19

"""

"""
}


process Macosko_pl_mapping {

publishDir params.outdir, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /${name}$/) "STAR_mapping/$filename"
}

input:
 set name, paramfile from g_7_paramset_g_19

output:
 file "${name}" into g_19_folder
 file "${name}" into g_19_folder2_g_20

"""
# Read paramfile into a var
paramline=\$(<${paramfile})

# Split line
IFS=\$'\t'
params=(\${paramline})

# Create param vars
read1=\${params[1]}
read2=\${params[2]}
#nestcells=\${params[3]}
genomedir=\${params[3]}
referencefasta=\${params[4]}

#echo \${read1} \${read2} \${genomedir} \${referencefasta} > ${name}.out.params.txt

# Define STAT exec path
STAR="/PATH/TO/STAR"

# add bc # Added tmp dir
cmd_pre="java -Djava.io.tmpdir=./tmpdir -jar /PATH/TO/Drop-seq_tools-2.0.0/3rdParty/picard/picard.jar FastqToSam \
F1=\${read1} \
F2=\${read2} \
SM=${name} \
O=${name}.unmapped.bam"

# Execute
eval \${cmd_pre}

# Create out dir
mkdir ${name}
# use -k to keep intermediates otherwise there is an error.

cmd_map="bash /PATH/TO/Drop-seq_tools-2.0.0/Drop-seq_alignment.sh \
-g \${genomedir} \
-r \${referencefasta} \
-o ${name} \
-d /PATH/TO/bin/Drop-seq_tools-2.0.0 \
-t ${name} \
-s \${STAR} -k \
./${name}.unmapped.bam"

# Execute
eval \${cmd_map}

# Write commands to file
echo \${cmd_pre} > ./${name}/exec.pre.cmd.txt
echo \${cmd_map} > ./${name}/exec.map.cmd.txt

# Write STAR version to file
\${STAR} | head -8 > ./${name}/STAR.version

# Copy the param file to the out directory
cp ${paramfile} ./${name}/params.txt

# Copy the unmapped bam file to out dir
cp ./${name}.unmapped.bam ./${name}

# final bam is now named ${name}/final.bam
# Rename to ${name}.bam
mv ./${name}/final.bam ./${name}/${name}.bam
mv ./${name}/final.bai ./${name}/${name}.bai

# RENAME the output file
#mv ./${name}/error_detected.bam ./${name}/${name}.bam
#mv ./${name}/error_detected.bai ./${name}/${name}.bai



"""
}


process DigitalExpression {

publishDir params.outdir, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /${name}$/) "count_mat/$filename"
}

input:
 file name from g_19_folder2_g_20

output:
 file "${name}" into g_20_folder
 set val("${name}"), file('mapdata') into g_20_paramset_g_22

"""

# Create folder for full matrix
mkdir mat_full
cmd="/PATH/TO/bin/Drop-seq_tools-2.0.0/DigitalExpression \
I=${name}/${name}.bam \
O=mat_full/${name}.dge.txt.gz \
SUMMARY=mat_full/${name}.dge.summary.txt \
TMP_DIR=mat_full \
MIN_NUM_GENES_PER_CELL=1"


# Execute
eval \${cmd}

# Write command to file
echo \${cmd} > mat_full/exec.cmd.txt

# Create reduced mat
mkdir mat_min_200_genes_per_cell
cmd2="/PATH/TO/bin/Drop-seq_tools-2.0.0/DigitalExpression \
I=${name}/${name}.bam \
O=mat_min_200_genes_per_cell/${name}.dge.txt.gz \
SUMMARY=mat_min_200_genes_per_cell/${name}.dge.summary.txt \
TMP_DIR=mat_min_200_genes_per_cell \
MIN_NUM_GENES_PER_CELL=200"
# Execute
eval \${cmd2}

# Write command to file
echo \${cmd2} > mat_min_200_genes_per_cell/exec.cmd.txt

# Move data dir to mapdata
mv ${name} mapdata

# Create results dir and mv results
mkdir ${name}
mv mat_full ${name}
mv mat_min_200_genes_per_cell ${name}

# Create the barcodes file
# cut -f1 ${name}/mat_min_200_genes_per_cell/${name}.dge.summary.txt | tail -n +3 > ${name}/mat_min_200_genes_per_cell/barcodes_min_200.txt
# Print barcodes (col 1) after matching "NUM_GENIC_READS"
sed -e '1,/NUM_GENIC_READS/d'  ${name}/mat_min_200_genes_per_cell/${name}.dge.summary.txt | cut -f1 > ${name}/mat_min_200_genes_per_cell/barcodes_min_200.txt

# Copy barcodesfile also to mapdata folder
cp ${name}/mat_min_200_genes_per_cell/barcodes_min_200.txt mapdata




"""
}


process Velocyto_folder_dropest {

publishDir params.outdir, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /${name}$/) "velocyto/$filename"
}

input:
 val ncpus from g_16_ncpus_g_22
 set name, mapdata from g_20_paramset_g_22

output:
 file "${name}" into g_22_folder

"""
export LC_ALL=C.UTF-8
export LANG=C.UTF-8

# Read the paramsfile and split to get gft file
# Read paramfile into a var
paramline=\$(<${mapdata}/params.txt)

# Split line
IFS=\$'\t'
params=(\${paramline})

# gtf is 6. param (note 0 based indexing)
gtffile=\${params[5]}

# Get version
velocyto --version > version.txt
vcversion=`egrep -o "([0-9]{1,}\\.)+[0-9]{1,}" version.txt`

# Copy the bc file to cwd
cp ${mapdata}/barcodes_min_200.txt .

# Create a link to bamfile
ln -s ${mapdata}/${name}.bam input.bam

## filter the bam file by the barcodes and add CB tag
## Based on script of Matthias Heinig
samtools view input.bam | \
awk -v barcode_file=barcodes_min_200.txt '
BEGIN{
  ## first read the barcodes that correspond to cells
  barcode = ""
  while ((getline line < barcode_file) > 0)
    barcode = barcode " " line
  close(barcode_file)
  nbarcodes = split(barcode, barcode_arr1, " ")
  ## the "in" function works on the hash keys so we need a second array
  ## with the barcodes as hash keys
  for (i in barcode_arr1) {
    barcode_arr2[barcode_arr1[i]] = barcode_arr1[i]
  }
  print "Considering " nbarcodes " barcodes" > "/dev/stderr"
}{
  # Barcode format XC:Z:ACCCTGAACACT
  match(\$0, /XC:Z:[ACGT]{12}/)
  if (RSTART > 0) { 
    barcode = substr(\$0, RSTART + 5, 12)
    if (barcode in barcode_arr2) {
      ## add the barcode as CB tag
      out = \$0 "\tCB:Z:" barcode
	  # out = \$0
      print out 
    }
  }
}' > aln_filtered.sam

# eval \${filtercmd}

## Get original header
samtools view -H ${mapdata}/${name}.bam > header

## create an indexed bam file
cat header aln_filtered.sam | samtools view -S -b -o aln_filtered.bam -
samtools index aln_filtered.bam

## remove tmp files
rm -f aln_with_readgroups.sam header

# Sort bamfile by CB tag / for velocyto
samtools sort --threads ${ncpus} -t CB -O BAM -o cellsorted_aln_filtered.bam aln_filtered.bam

# Run velocyto
# velocmd="velocyto run --samtools-threads ${ncpus} -b ${mapdata}/barcodes_min_200.txt -o velocyto ${mapdata}/${name}.bam \${gtffile}"
# Run velocyto on filtered bam file
# Velocyto will look for the file cellsorted_aln_filtered.bam
velocmd="velocyto run --samtools-threads ${ncpus} -o velocyto aln_filtered.bam \${gtffile}"
eval \${velocmd}

# mv the original folder (which is a softlink)
# mv ${name} oldname

# mv velocyto folder to name
mv velocyto ${name}

# cp velocyto version file
cp version.txt ${name}/vc_version.txt
# Save velocyto command
echo \${velocmd} > ${name}/vc.cmd.txt



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
