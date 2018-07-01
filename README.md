# MDIBL-T3-WGS
Bacterial Genome Assembly and Assessment Tutorial


## Various Resources:
[MDIBL T3 Course Website](https://labcentral.mdibl.org/workspaces/view/5ad10ee2-cf8c-4894-a980-1172d1dec312/pages/5ad8c98a-76a8-4b42-a40d-18a4d1dec312)

[Beginning Bash Cheatsheet](https://maker.pro/linux/tutorial/basic-linux-commands-for-beginners)

[Filezilla Download](https://filezilla-project.org/download.php)

[PuTTY Download](https://www.putty.org/)


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
#
```

* Examine Raw Reads
       - Note the extension, this data is compressed.
```bash
# 

```
![rawilluminadatafastqfiles](https://user-images.githubusercontent.com/18738632/42129269-49b8dace-7c8e-11e8-86e7-069df9028447.png)

* Count The Number of Raw Reads
```bash
# 
```

* Whats our total bp of data?
* If we have a 5 MB genome, what is our average coverage?

## Read Quality Check w/ FASTQC
manual: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

alterative tools: Just use fastqc

[FASTQC explained](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/)

* Run Fastqc
```bash
#
```

* Transfer resulting HTML files to computer using filezilla or with the commandline on OSX/Linux.
```bash
# In a fresh terminal on OSX or Linux
scp USERNAME@ron.sr.unh.edu:/home/maineBK/USERNAME/mdibl-t3-2018-WGS/fastqc_raw-reads/*.html /path/to/put/files
```

## Adapter and Quality Trimming w/ Trimmomatic
manual: http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf

alternative tools: [cutadapt](http://cutadapt.readthedocs.io/en/stable/guide.html),[skewer](https://github.com/relipmoc/skewer)

* Run Trimmomatic
```bash
# 
```
* Move trimmed reads to new directory
```bash
#
```

* How do I know where scripts are held?
```bash
# 
```
* Run FastQC again.
* Count the number of reads in the new files

## Genome Assembly w/ SPAdes
manual: http://cab.spbu.ru/software/spades/

alternative tools: [ABySS](http://www.bcgsc.ca/platform/bioinfo/software/abyss),[MaSuRCA](http://masurca.blogspot.com/)

* Run SPAdes
```bash
#
```
* View Output Data
```bash
# 
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

alternative tools: [NCBI PGA](https://www.ncbi.nlm.nih.gov/genbank/genomesubmit_annotation/)

* Run PROKKA
```bash
#
```
![prokka_workflow](https://user-images.githubusercontent.com/18738632/42130490-e45251b6-7cb4-11e8-99ef-9579b9b7ce05.png)

![gene_annotatoion](https://user-images.githubusercontent.com/18738632/42130642-bf1fb57e-7cb8-11e8-8472-37b82dadb53e.png)

