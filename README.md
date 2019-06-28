# MDIBL-T3-WGS
Bacterial Genome Assembly and Assessment Tutorial

**This site will be continuously updated throughout the week.**

## General Overview

   Throughout this tutorial we will be going through the process of *de novo* genome assembly. This workflow begins with raw sequencing data you would recieve from a standard sequencing center (in fastq format). We start by examining the fastqs for quality with **fastqc**. Next we trim the low quality bases and remove adapter sequences from the reads with **Trimmomatic**. Once that is done we move directly into genome assembly with **SPAdes**. The SPAdes pipeline does the brunt of the work, taking our trimmed sequencing reads as input and providing a FASTA file as output, this FASTA file is our genome assembly. From here we assess the genome assembly for contiguity using **QUAST** and for content/comleteness with **BUSCO**. Finally we use several different programs including **BLAST**, **BWA**, and **blobtools**, to filter the genome for potential contaminates/non-target sequences. At this point you should have a novel genome that is ready for submission to NCBI and/or for comparative genomics with previously published genomes.

## Various Resources:
[MDIBL T3 Course Website](https://labcentral.mdibl.org/workspaces/view/5ad10ee2-cf8c-4894-a980-1172d1dec312/pages/5ad8c98a-76a8-4b42-a40d-18a4d1dec312)

[Beginning Bash Cheatsheet](https://maker.pro/linux/tutorial/basic-linux-commands-for-beginners)

[BASH tutorials](http://nhinbre.org/wp-content/uploads/2015/06/intro_worksheet-1.pdf)

[Filezilla Download](https://filezilla-project.org/download.php)

[PuTTY Download](https://www.putty.org/)

Any Questions when working through this? Email me: jlsevigny1@wildcats.unh.edu
