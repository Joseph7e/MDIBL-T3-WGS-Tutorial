# MDIBL-T3-WGS
Bacterial Genome Assembly and Assessment Tutorial

![alt text](https://www.digitalcheck.com/wp-content/uploads/2018/08/old-compouter-zenith-terminal-400px.jpg)

**This site will be continuously updated throughout the week.**

### General Notes:
**For each program that we run in this tutorial I have provided a link to the manual**. These manuals provide a thorough explanation of what exactly we are doing. Before running the program it is a good idea to skim through these, examine the options, and see what it does. It is also a good idea to check out the publication associated with the program. Please note that the commands we run are general and usually executed with default settings. This works great for most genomes but the options may need to be tweaked depending on your genome. Before you run any command it is also a great idea to look at the programs help menu. This can usually be done with the name of the program followed by '-h' or '-help' or '--help'. i.e. 'spades -h'. Also ... never forget about google for quick answers to any confusion.

This tutorial assumes a general understanding of the BASH environment. **You should be familiar with moving around the directories and understand how to manipulate files**.

Throughout this tutorial the commands you will type are formatted into the gray text boxes (don't do it when learning but they can be faithfully copied and pasted). The '#' symbol indicates a comment, BASH knows to ignore these lines. 

**Remember to tab complete!** There is a reason the tab is my favorite key. It prevents spelling errors and allows you to work 10X faster (I timed it). Remember if a filename isn't auto-completing you can hit tab twice to see your files while you continue typing your command. If a file doesn't auto-complete it means you either have a spelling mistake, are in a different directory than you originally thought, or that it doesn't exist.

## Starting Data:
Your starting data is found within a shared directory within your group folder (one directory level up). To start we will move a set of Sample data into your home directories. Each of these samples represent the genome of a unique and novel microbe that has not been seen before (except by me). Inside this directory are Illumina HiSeq 2500, paired-end, 250 bp sequencing reads. Looking in this directory you should see two files per sample, the forward and reverse reads. These files are in **FASTQ** format (see below).

* Get your bearing on the server.

It's hard to know where your going if you don't know where you are. When I am working on the server I constantly type 'ls' and 'pwd' to make sure I am where I think I am. You should too!

```bash
# print your current working directory. If you just logged in you should be in your home directory (/home/group/username/)
pwd
# change to your home directory in case you weren't already there. Remember ~/ is an absolute path to your home directory.
cd ~/
# ls to view your home dir
ls
# ls to view the group directory (I type ‘ls’ a lot!)
ls ../
# view the shared directory of starting data
ls ../shared
# View the shared project with the ‘tree’ command
tree ../shared
# Copy a sample from the shared directory to your home dir, 
#“Project_X”, where X denotes the Project name and "Sample_X" (where X denotes your sample name).
# USE AUTOCOMPLETE
cp -r ../shared/Project_X/ Sample_X/ ./
# confirm the copy arrived (remember ‘*’ will match any character/string)
ls Sample_*/
```

[Link explaining the 'Read Name Format'](http://support.illumina.com/content/dam/illumina-support/help/BaseSpaceHelp_v2/Content/Vault/Informatics/Sequencing_Analysis/BS/swSEQ_mBS_FASTQFiles.htm): SampleName_Barcode_LaneNumber_001.fastq.gz


Important note: In the above command I use the "\*" character to view the Sample directory, I would normally just type out the entire path using tab complete (which is what you should do). This wildcard will match any string of characters. I use this because everyone will have a different Sample name. To make this tutorial as general as possible I need to use these wildcards throughout the tutorial. In addition I may use Sample_X instead of Sample_\*. In these cases be sure to type out your complete sample name!, the wildcards probably won't work 


* Prepare your working directory

It is a good idea to keep your directories tidy and to name your files something that makes sence. This is just to keep things organized so you know what everything is several months from now. We are going to make a new directory to house all of the analyses for this tutorial.

```bash
# Make a new directory and add the Sample directory into it
mkdir mdibl-t3-2019-WGS
mv Sample* mdibl-t3-20189-WGS/
cd mdibl-t3-2019-WGS/
# make the sample directory name more meaningful
mv Sample_X Sample_X-raw_reads
```

## Examine the Raw Reads

Note the file extension - fastq.**gz**. Since these files are usually pretty big it is standard to receive them compressed. To view these files ourselves (which you normally wouldn't do) you either have to decompress the data with gzip or by using variations of the typical commands. Instead of 'cat' we use 'zcat', instead of grep we can use 'zgrep'. Below I show both ways.
       
```bash
# Examine the reads with zcat, I use the wildcard '*' to match the file since everyone's names will be different. Use tab complete and there is no need for the wildcards.
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

    - Line 1. Always begins with an '@' symbol and denotes the header. This is unique to each sequence and has info about the sequencing run. 

    - Line 2. The actual sequencing read for your organism, a 250 bp string of As, Ts, Cs, and Gs.

    - Line 3. Begins with a '+' symbol, this is the header for the read quality. Usually the same as the first line header. 

    - Line 4. Next are ascii symbols representing the quality score (see table below) for each base in your sequence. This denotes how confident we are in the base call for each respective nucleotide. This line is the same length as the sequencing line since we have a quality score for each and every base of the sequence. 

![rawilluminadatafastqfiles](https://user-images.githubusercontent.com/18738632/42129269-49b8dace-7c8e-11e8-86e7-069df9028447.png)

![quality_info](https://user-images.githubusercontent.com/18738632/42226531-2f343178-7ead-11e8-8401-5a2fb455b4ef.png)



* Count The Number of Raw Reads

I always start by counting the number of reads I have for each sample. This is done to quickly assess whether we have enough data to assemble a meaningful genome. Usually these file contains millions of reads, good thing BASH is great for parsing large files! Note that the forward and reverse reads will have the same number of entries so you only need to count one.

```bash
# using grep. Note that I don't count just '@', this is because that symbol may appear in the quality lines.
zgrep -c '@HSQ' Sample*/*R1*
# counting the lines and dividing by 4. Remember each read entry is exactly four lines long. These numbers should match.
zcat Sample*/*_R1_* | wc -l
```
* Whats our total bp of data? This is what we call our sequencing throughput. We multiple the number of reads by the read length (ours is 250) and by 2 because it is paired-end data.

(Read length x 2(paired-end) x Number of reads)

```
# we can do this calculation from the terminal with echo and bc (bc is the terminal calculator)
echo "Number_of_reads * 250 * 2" | bc
```

* If we have a 7 MB genome, what is our average coverage? 

(Total bp/7,000,000)
```
echo "Total_bps / 7000000" | bc
```

If you completed the above calculation lets hope you have at least 10X coverage. For the most part, the higher the coverage the better off we are. If you have low coverage you'll want to do some more sequencing and get more read data. Usually published genomes have at least 70-100X coverage.

## Read Quality Check w/ FASTQC
manual: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

[FASTQC explained](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/)

* Run Fastqc

FastQC is a program to summarize read qualities and base composition. Since we have millions of reads there is no practical way to do this by hand. We call the program to parse through the fastq files and do the hard work for us. **The input to the program is one or more fastq file(s) and the output is an html file with several figures.** The link above describes what each of the output figures are describing. I mainly focus on the first graph which visualizes our average read qualities and the last figure which shows the adapter content. Note that this program does not do anything to your data, as with the majority of the assessment tools, it merely reads it.

```bash
# make a directory to store the output
mkdir fastqc_raw-reads
# run the program
fastqc Sample_*/*_R1_* Sample_*/*_R2_* -o fastqc_raw-reads
ls fastqc_raw-reads
# the resulting folder should contain a zipped archive and an html file, we can ignore the zipped archive which is redundant.
```

* Transfer resulting HTML files to computer using filezilla or with the command line on OSX/Linux.

On filezilla you will need to enter the same server information when you login form the terminal. Be sure to use port 22.  

```bash
# to get the absolute path to a file you can use the ‘readlink’ command.
readlink -f fastqc_raw-reads/*.html
# copy those paths, we will use them for the file transfer
# In a fresh terminal on OSX, Linux, or BASH for windows
scp USERNAME@ron.sr.unh.edu:/home/GROUP/USERNAME/mdibl-t3-2019-WGS/fastqc_raw-reads/*.html /path/to/put/files
```
