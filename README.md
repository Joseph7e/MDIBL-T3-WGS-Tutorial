# MDIBL-T3-WGS
Bacterial Genome Assembly and Assessment Tutorial

       Throughout this tutorial we are going through the process of *de novo* genome assembly. This process begins with raw sequencing data in the form of fastqs you would recieve from a sequencing center. We start by examing the fastqs for quality with **fastqc**. Next we trim low quality bases from our reads and remove adapter sequences with **Trimmomatic**. Once we are happy with the quality of the bases we move directly into genome assembly with **SPAdes**. This program takes our trimmed sequencing reads and provides a FASTA file, this is our genome assembly. From here we assess the genome assembly for contiguity using **QUAST** and for content/comleteness with **BUSCO**. From there we use several different programs including **BLAST**, **BWA**, and **blobtools**, to filter the genome for potential contaminates/non-target sequences. 




## Various Resources:
[MDIBL T3 Course Website](https://labcentral.mdibl.org/workspaces/view/5ad10ee2-cf8c-4894-a980-1172d1dec312/pages/5ad8c98a-76a8-4b42-a40d-18a4d1dec312)

[Beginning Bash Cheatsheet](https://maker.pro/linux/tutorial/basic-linux-commands-for-beginners)

[Filezilla Download](https://filezilla-project.org/download.php)

[PuTTY Download](https://www.putty.org/)


### General Notes:
For each program that we run there are links to the manuals, some more helpful than others. These manuals provide a thorough explanation of what exactly we are doing so be sure to skim through them. 

## Starting Data:
Your starting data is in a directory called "Sample_X" (where X donates your sample name). Inside this directory are Illumina HiSeq 2500, paired-end, 250 bp sequencing reads.

* Get your bearing on the server and view data
```bash
# print your current working directory
pwd
# ls to view your read directory
ls Sample*
# use the 'tree' command
tree
tree -L 2
```
[Read Name Format](http://support.illumina.com/content/dam/illumina-support/help/BaseSpaceHelp_v2/Content/Vault/Informatics/Sequencing_Analysis/BS/swSEQ_mBS_FASTQFiles.htm): SampleName_Barcode_LaneNumber_001.fastq.gz


* Prepare working directory - keep things organized for your future self
```bash
# Make a new directory and add Sample directory into it
mkdir mdibl-t3-2018-WGS
mv Sample* mdibl-t3-2018-WGS/
cd mdibl-t3-2018-WGS/
# make the sample directory name more meaningful
mv Sample_X Sample_X-raw_reads
```

* Examine Raw Reads
       - Note the extension, this data is compressed.
```bash
# Examine the reads with zcat
zcat Sample*/*_R1_* | more
# unzip the data and view
gunzip Sample*/*_R1_*
more Sample*/*_R1_*
# rezip the data
gzip Sample*/*_R1_*
# Examine reads directly with less
less -S Sample*/*_R1_*
```
![rawilluminadatafastqfiles](https://user-images.githubusercontent.com/18738632/42129269-49b8dace-7c8e-11e8-86e7-069df9028447.png)

* Count The Number of Raw Reads
```bash
# using grep
zgrep -c '@HSQ' Sample*/*R1*
# counting the lines and divide by 4
zcat Sample*/*_R1_* | wc -l
```

* Whats our total bp of data? (Read length x 2(paired-end) x Number of reads)
* If we have a 5 MB genome, what is our average coverage? (Total bp/5,000,000)

## Read Quality Check w/ FASTQC
manual: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

alterative tools: Just use fastqc

[FASTQC explained](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/)

* Run Fastqc
```bash
mkdir fastqc_raw-reads
fastqc Sample_*/*_R1_* Sample_*/*_R2_* -o fastqc_raw-reads
ls fastqc_raw-reads
# the resulting folder contains a zipped archive and an html file
```

* Transfer resulting HTML files to computer using filezilla or with the commandline on OSX/Linux.
```bash
# In a fresh terminal on OSX or Linux
scp USERNAME@ron.sr.unh.edu:/home/maineBK/USERNAME/mdibl-t3-2018-WGS/fastqc_raw-reads/*.html /path/to/put/files
```

## Adapter and Quality Trimming w/ Trimmomatic
manual: http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf

alternative tools: [cutadapt](http://cutadapt.readthedocs.io/en/stable/guide.html), [skewer](https://github.com/relipmoc/skewer)

* Run Trimmomatic
```bash
# Run wrapper script
trim_script_TruSeq.sh <forward_reads> <reverse_reads>
```
* Move trimmed reads to new directory
```bash
mkdir trimmed_reads
mv paired* unpaired* trimmed_reads/ 
```

* How do I know where scripts are held?
```bash
# determine where the script is held
which trim_script_TruSeq.sh
# view it
more /usr/local/bin/trim_script_TruSeq.sh
# examine adapters, path found in above command
more /usr/share/trimmomatic/TruSeq2-PE.fa
```
* Run FastQC again.
* Count the number of reads in the new files (see above)

## Genome Assembly w/ SPAdes
manual: http://cab.spbu.ru/software/spades/

alternative tools: [ABySS](http://www.bcgsc.ca/platform/bioinfo/software/abyss),[MaSuRCA](http://masurca.blogspot.com/)

* Run SPAdes
```bash
# examine the help menu
spades.py --help
# using "nohup some_command &" allows you to run the job on the server while your laptop is closed/off 
nohup spades.py -1 trimmed_reads/paired_forward.fastq.gz -2 trimmed_reads/paired_reverse.fastq.gz -s trimmed_reads/unpaired_forward.fastq.gz -s trimmed_reads/unpaired_reverse.fastq.gz -o spades_assembly_default &
```
* View Output Data
```bash
ls spades_assembly_default/
# view FASTA file
less -S spades_assembly_default/contigs.fasta
# view the top 10 headers
grep '>' spades_assembly_default/contigs.fasta | head
```

* Clean up Spades dictory.
It is a good idea to keep your directories clean and organized. A part of this means cleaning up unwanted files.
    
```bash
# 
```

* How can I remove all the unwanted data and keep the contigs and log?
* How can I view just the headers? How about the top 10?
* How can I count the number of contigs I have?

## Genome Structure Assessment w/ QUAST
manual: http://quast.bioinf.spbau.ru/manual.html

* Run Quast
```bash
#
```
* View ouput files
```bash
#
```

## Genome Content Assessment w/ BUSCO
manual: https://busco.ezlab.org/

* BUSCO preperation
```bash
# Busco requires a path variable to be set before use.
echo export AUGUSTUS_CONFIG_PATH="/usr/local/src/augustus-3.2.2/config/" >> ~/.bashrc
source ~/.bashrc
```

* Path to lineage data on RON: 
```bash
# View available sets
ls /usr/local/src/augustus-3.2.2/rc/
```
* Run BUSCO
```bash
#
```

## Genome Annotation w/ PROKKA
manual: https://github.com/tseemann/prokka

alternative tools: [NCBI PGA](https://www.ncbi.nlm.nih.gov/genbank/genomesubmit_annotation/), [Glimmer](https://ccb.jhu.edu/software/glimmer/)

* Run PROKKA
```bash
grep -o "product=.*" prokka_output/PROKKA_07022018.gff | sed 's/product=//g' | sort | uniq -c | sort -nr > protein_abundances.txt
```



* Extract the 16S sequence from the FFN file.



![prokka_workflow](https://user-images.githubusercontent.com/18738632/42130490-e45251b6-7cb4-11e8-99ef-9579b9b7ce05.png)

![gene_annotatoion](https://user-images.githubusercontent.com/18738632/42130642-bf1fb57e-7cb8-11e8-8472-37b82dadb53e.png)



## NCBI BLAST


## Command Line BLAST

* Make a BLAST db from your contig files


## Read Mapping w/ BWA and samtools

* Full workflow
```bash
bwa index -a bwtsw $fasta
bwa mem -M -t 24 $fasta $forward $reverse > raw_mapped.sam
samtools view -@ 24 -Sb -F 4  raw_mapped.sam  | samtools sort -@ 24 - -o sorted_mapped.bam
#bioawk -c fastx '{ print $name, length($seq) }' < $fasta > genome.txt
#genomeCoverageBed -ibam sorted_mapped.bam -g genome.txt > coverage.txt
bedtools genomecov -ibam sorted_mapped.bam > coverage.out
gen_input_table.py  --isbedfiles $fasta coverage.out >  coverage_table.tsv
```

## Non-target contig removal w/ Blobtools

* Generate "Hits" file

* RUnnning BLAST on the server
```bash
blobtools create -i contigs.fasta -b raw_mapped.sam -t contigs.fasta.vs.nt.cul5.1e5.megablast.out -o blob_out
blobtools view -i blob_out.blobDB.json -r all -o blob_taxonomy
blobtools plot -i blob_out.blobDB.json -r genus

```
