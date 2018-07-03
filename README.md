# MDIBL-T3-WGS
Bacterial Genome Assembly and Assessment Tutorial

## General Overview

   Throughout this tutorial we will be going through the process of *de novo* genome assembly. This workflow begins with raw sequencing data, exactly how you would recieve it from a standard sequencing center (fastqs). We start by examining the fastqs for quality with **fastqc**. Next we trim the low quality bases and remove adapter sequences from the reads with **Trimmomatic**. Once that is done we move directly into genome assembly with **SPAdes**. This program does the brunt of the work, taking our trimmed sequencing reads as input and providing a FASTA file, this FASTA file is our genome assembly. From here we assess the genome assembly for contiguity using **QUAST** and for content/comleteness with **BUSCO**. Finally we use several different programs including **BLAST**, **BWA**, and **blobtools**, to filter the genome for potential contaminates/non-target sequences. At this point you should have a novel genome that is ready for submission to NCBI and for comparative genomics with previously published genomes.
   

## Various Resources:
[MDIBL T3 Course Website](https://labcentral.mdibl.org/workspaces/view/5ad10ee2-cf8c-4894-a980-1172d1dec312/pages/5ad8c98a-76a8-4b42-a40d-18a4d1dec312)

[Beginning Bash Cheatsheet](https://maker.pro/linux/tutorial/basic-linux-commands-for-beginners)

[More BASH tutorials] (http://nhinbre.org/wp-content/uploads/2015/06/intro_worksheet-1.pdf)

[Filezilla Download](https://filezilla-project.org/download.php)

[PuTTY Download](https://www.putty.org/)


### General Notes:
For each program that we run there are links to the manuals. These manuals provide a thorough explanation of what exactly we are doing. It is important to at least skim through these to examine the options and what it does. The commands we run are usually general and rely on default settings, this works great for most genomes but the options may need to be tweaked for each genome. Before you run any command it is also a great idea to look at the programs help menu. This can usually be done with the name of the program followed by '-h' or '-help' or '--help'. i.e. **spades -h**. Also ... never forget about google for quick answers to any confusion.

Also note that this tutorial assumes a general understanding of the BASH environment. You should be familiar with moving around the directories and understand how to manipulate files.

Throughout the tutorial the commands you will type are formatted into the grey text boxes and can be copied and pasted. The '#' symbol indicates a comment and BASH knows to ignore these lines. 

**Remember to tab complete!** There is a reason the tab is my favorite key. It prevents spelling errors and allows you to work 10X faster (I timed it). Remember if a filename isn't auto-completing you can hit tab twice to see your file options while you continue typing your command. If the file still doesn't auto-complete it means you either have a spelling mistake or are in a different directory than you originally thought.

## Starting Data:
Your starting data is in a directory called "Sample_X" (where X donates your sample name). I placed a different Sample directory into each of your directories, each represents a unique and novel microbe that has not been seen before (except by me). Inside this directory are Illumina HiSeq 2500, paired-end, 250 bp sequencing reads. Looking in this directory you should see two files per sample, the forward and reverse reads. These files are in **FASTQ** format (see below). 

* Get your bearing on the server - it's hard to know where your going if you don't know where you are.
```bash
# print your current working directory. If you just logged in you should be in your home directory (/home/group/username/)
pwd
# change to your home directory in case you weren't already there. Remember ~/ is an absolute path to your home directory.
cd ~/
# ls to view your read directory.
ls Sample_*
# use the 'tree' command to view your current directory structure.
tree
# if things are messy you can limit the number of sub-directories you see with the '-L' option.
tree -L 2
```

** Quick note: In the above command I use the "\*" charcter to view the Sample directory. This is a wildcard which means match any character. I use this because everyone will have a different Sample name. In commands below I may use Sample_X instead of Sample_\*. In these cases be sure to type out your complete sample name!  ** 

[Link explaining the 'Read Name Format'](http://support.illumina.com/content/dam/illumina-support/help/BaseSpaceHelp_v2/Content/Vault/Informatics/Sequencing_Analysis/BS/swSEQ_mBS_FASTQFiles.htm): SampleName_Barcode_LaneNumber_001.fastq.gz


* Prepare your working directory - It is a good idea to keep your directories tidy and to name your files somehting that makes since. This is just to keep things organized so it makes sense to your future self.
```bash
# Make a new directory and add Sample directory into it
mkdir mdibl-t3-2018-WGS
mv Sample* mdibl-t3-2018-WGS/
cd mdibl-t3-2018-WGS/
# make the sample directory name more meaningful
mv Sample_X Sample_X-raw_reads
```

* Examine the Raw Reads

Note the file extension - fastq.**gz**. Since these files are usually pretty big it is standard to recieve them compressed. To view these files ourselves (which you normally wouldn't do) you either have to decompress the data with gzip or by using variations of the typical commands. Instead of 'cat' we use 'zcat', instead of grep we can use 'zgrep'. Below I show both ways.
       
```bash
# Examine the reads with zcat, I use the wildcard '*' to match the file since everyones names will be different. Use tab complete and there is no need for the wildcards.
zcat Sample*/*_R1_* | more
# unzip the data and view
gunzip Sample*/*_R1_*
more Sample*/*_R1_*
# rezip the data
gzip Sample*/*_R1_*
# Examine reads directly with less
less -S Sample*/*_R1_*
```

* Fastq File Format - each sequencing read entry is four lines long.. 

    - Line 1. Always begins with an '@' symbol and donates the header. This is unique to each sequence and has info about the sequncing run. 

    - Line 2. The actual sequencing read for your organism, a 250 bp string of As, Ts, Cs, and Gs.

    - Line 3. Begins with a '+' symbol, this is the header for the read quality. Usually the same as the first line header. 

    - Line 4. Next are ascii symbols representing the quality score (see table below) for each base in your sequence. TThis donates how confident we are in the base call for each respective nucleotide. This line is the same length as the sequencing line since we have a quality score for each and every base of the sequence. 

![rawilluminadatafastqfiles](https://user-images.githubusercontent.com/18738632/42129269-49b8dace-7c8e-11e8-86e7-069df9028447.png)

![quality_info](https://user-images.githubusercontent.com/18738632/42226531-2f343178-7ead-11e8-8401-5a2fb455b4ef.png)



* Count The Number of Raw Reads

I always start by counting the number of reads I have for each sample. This is done to make sure we have enough data to assemble a meaningful genome in the first place. Usually the file contains millions of reads, good thing BASH is great for parsing large data files! Note that the forward and reverse reads will have the same number of entries so you only need to count one (at least they should!)

```bash
# using grep
zgrep -c '@HSQ' Sample*/*R1*
# counting the lines and divide by 4
zcat Sample*/*_R1_* | wc -l
```
* Whats our total bp of data? This is what we call our sequencing throughput 

(Read length x 2(paired-end) x Number of reads)

* If we have a 7 MB genome, what is our average coverage? 

(Total bp/7,000,000)

If you completed the above calculation lets hope you have at least 10X coverage. For the most part, the higher the coverage the better off we are. If you have low coverage you'll want to do some more sequencing and get more read data. Usually published genomes have at least 70-100X coverage.

## Read Quality Check w/ FASTQC
manual: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

alterative tools: Just use fastqc

[FASTQC explained](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/)

* Run Fastqc

Fastqc is a program to summarize the read qualities. Since we have millions of reads there is no practical way to do this by hand. We call the program to parse through the fastq files and do the hard work for us. **The input to the program is one or more fastq file(s) and the output is an html file with several figures.** The link above describes what each of the figures are showing. I mainly use the first figure which is our read qualities and the last figure which shows what sort of adapter content we have. Note that this program does not do anything to your data, it mearly reads it.

```bash
# make a directory to store the output
mkdir fastqc_raw-reads
# run the program
fastqc Sample_*/*_R1_* Sample_*/*_R2_* -o fastqc_raw-reads
ls fastqc_raw-reads
# the resulting folder should contain a zipped archive and an html file
```

* Transfer resulting HTML files to computer using filezilla or with the commandline on OSX/Linux.

On filezilla you will need to enter the same server information when you login form the terminal. Be sure to use port 22.  

```bash
# In a fresh terminal on OSX or Linux
scp USERNAME@ron.sr.unh.edu:/home/maineBK/USERNAME/mdibl-t3-2018-WGS/fastqc_raw-reads/*.html /path/to/put/files
```

## Adapter and Quality Trimming w/ Trimmomatic
manual: http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf

alternative tools: [cutadapt](http://cutadapt.readthedocs.io/en/stable/guide.html), [skewer](https://github.com/relipmoc/skewer)

* Run Trimmomatic

You may have noticed from the fastqc output the some of your reads had poor qualities towards the end of the sequences, this is especially true for the reverse reads. You may also notice that fastqc failed for adapter content. This programs will be used to trim these low quality bases and to remove the adapters. I created a wrapper script called trim_script_TruSeq.sh which makes this progam much easier to use. It is available on the server by calling its name, it is also available on this github repository. For this wrapper script **the input is the raw forward and reverse reads and the output will be new trimmed fastq files** which we will use for genome assembly. When you are more comfortable using BASH you can call trimmomatic directly by using the manual or by copying the code from the provided script.

```bash
# Run wrapper script
trim_script_TruSeq.sh Sample_*/*_R1_* Sample_*/*_R2_*
```
* Move the new trimmed reads to new directory - remember its a good idea to keep the directory clean.
```bash
mkdir trimmed_reads
# move all the 
mv *fastq.gz trimmed_reads/
# confirm that the files have moved
ls trimmed_reads/
```

When the program finishes it outputs four files. paired_forward.fastq.gz, paired_reverse.fastq.gz, and two unpaired reads. These output files are cleaned reads which hopefully have only highly confident sequences and have no adapters. Some sequences will be lost entirely, some will lose a few bases off the ends, and some won't be trimmed at all. When a reverse read is lost but a forward read is maintained, the forward read will be written to the unpaired_forward.fastq.gz file (and vise-versa)

Similiar to above you can run FASTQC again with your new trimmed reads. Comparing the original html and the new one you should note the differences (higher quality and no adapters).

You can also count the number of reads for each of your files. How does this compare to the original count? What percentage of your reads did you lose?

## Genome Assembly w/ SPAdes
manual: http://cab.spbu.ru/software/spades/

alternative tools: [ABySS](http://www.bcgsc.ca/platform/bioinfo/software/abyss),[MaSuRCA](http://masurca.blogspot.com/)


There are many programs that are used for genome assembly. For the most part they are run the same. The input will be a set of sequencing reads in fastq format and the output will be a FASTA file which is the genome assembly.


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
