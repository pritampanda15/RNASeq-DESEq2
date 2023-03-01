#########################################
# Identifying differential gene expression and biased allele expression from raw unfiltered unaligned data
#		By Kimberly C. Olney 
#
# Last updated: August 30th, 2016
#
#--------------------------------------
# Overview:
# 			Identifying patterns of biased allele expression and differential gene expression 
#
# Samples information:
#			Update to fit your project 
#
#
# RNAseq data processing overview: 
#			Raw reads were mapped using STAR read aligner. This was done to achieve high sensitivity to both SNPs and, indels. Specifically, the STAR 2-pass method (Pär G Engström et al.) using the suggested protocol with the default parameters.  
#			The STAR 2-pass approach detects splice junctions in a first alignment run and then those are used to guide the final alignment. STAR uses genome index files that must be saved in unique directories.
#			The first pass alignment job is then executed and splice junctions are indentified. A new index is then created using splice junction information contained in the file SJ.out.tab from the first pass. 
#			The resulting index is then used to produce the final alignments. The above step produces a BAM file, which we then put through the suggested GATK's best practices steps: adding read group information, sorting, marking duplicates and indexing.
#			Then, using GATK tool SplitNCigarReads, we can split the reads into exon segments (getting rid of Ns but maintaining grouping information) and hard-clip any sequences overhanging into the intronic regions. 
#			Overhang regions that have mismatch greater than the threshold are clipped to reduce the number of the called false variants. AT SplitNCigarReads, mapping qualities are also reassigned because STAR assigns good alignments a MAPQ of 255 which is meaningless to GATK. 
#     To correct for this, GATK’s ReassignOneMappingQuality read filter, is used to reassign all good alignments to the default value of 60. To identify differential gene expression cuffdiff is run on all pairwise comparisons. 
#     To identify allele expression we will need to generate a variant call fille from the individual bam files. Variant calling is then accomplished using GATK's HaplotypeCaller tool. 
#     The HaplotypeCaller tool takes into account the information about intron-exon split regions that is embedded in the BAM file by SplitNCigarReads. The HaplotypeCaller tool will perform “dangling head merging” operations 
#     and avoids using soft-clipped bases to minimize false positive and false negative calls. From here we will filter the VCF using bcftools to remove indels and to only include biallelic sites. Next we will use
#     GATK's SelectVariants to slipt our merged VCF to have a VCF per individual. This is down to speed up down stream analysis as programs won't have to filter through a very large VCF containing all the individuals as this is not necessary for our needs. 
#     Using SnpEff's SnpSift we will filter the individual VCF to only include heterozygous sites as homozygous sites will be uninformative for inferring biased allele expression. Finally we will run GATK's ASEReadCounter to get allele count information for 
#     for each individual at each heterozygous site for the whole genome. 
#
#
#--------------------------------------
# Contents:
#			1. Download or obtain data 
#			2. convert data to be in fastq format
#			3. FastQC to check the quality of the raw reads
#			4. Trim fastq files for quality and to remove adaptors. 
#			5. FastQC to check the quality of the trimmed reads
#     6. Obtain reference genome and gene annotation files
#     7. Generate genome indexes
#			8. STAR: map transcript reads to the reference genome and identify splice junctions 
#     9. check quality of raw BAM files
#			10. Sort BAM files, to be in the same order as the reference for downstream analysis 
#     11. check quality of sorted BAM files
#			12. Mark duplicates, duplicates will not be removed but will be marked for quality checks
#     13. check quality of mark duplicates BAM files
#			14. Add read groups to BAM files, this done to keep the sample ids organized when creating the merged vcf file
#     15. check quality of add read group BAM files
#     16. Index BAM files
#			17. SplitNCigarReads, hard-clip any sequences overhanging into the intronic regions
#     18. Identify genes that differentially expressed 
# 	  19. Call variant calling using GATK's HaplotypeCaller tool 
# 		20. merge gvcf files to a raw group VCF file 
#     21. Filter the VCF to remove indels and to only keep heterozygous biallelic sites. 
#     22. Create 1 VCF per sample using SelectVariants
#     23. Filter for only heterozygous sites
#     24. Obtain allele count information using ASEReadCounter
#
#--------------------------------------
# Publicly available packages:
#       SRA toolkit                         https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?cmd=show&f=software&m=software&s=software
#       fastqc                              http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
#       Trimmomatic                         http://www.usadellab.org/cms/?page=trimmomatic
#		    STAR											          https://github.com/alexdobin/STAR
#       bamtools                            https://github.com/pezmaster31/bamtools
#		    Picard											        http://sourceforge.net/projects/picard/
#	     	GATK											          https://www.broadinstitute.org/gatk/download/auth?package=GATK
#       cuffdiff                            http://cole-trapnell-lab.github.io/cufflinks/cuffdiff/
#       bcftools                            https://samtools.github.io/bcftools/
#       SnpEff                              http://snpeff.sourceforge.net/download.html
#   
########################################
#--------------------------------------
# Install packages 
#-------------------------------------
# Packages installed on my home directory on NGCC in my tools folder  /tools
# GATK													        https://www.broadinstitute.org/gatk/download/auth?package=GATK
# STAR-master										        https://github.com/alexdobin/STAR
# Picard 												        http://sourceforge.net/projects/picard/
# Trimmomatic                           http://www.usadellab.org/cms/?page=trimmomatic
# SnpEff                                http://snpeff.sourceforge.net/download.html

########################################
#--------------------------------------
# 1. Download data
#--------------------------------------
# In sra (sequence read archive, as know as short-read archive) format 
# Include GEO accession number 
# Samples downloaded were saved in Project directory /Project/
# 
# Example command: $wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP/SRR1850937.sra
# wget                                                          stands for "web get"
# ftp                                                           stands for File transfer program
# sampleID.sra                                                  path to where the samples are located and the sampleID

# SBATCH script - 
sbatch NameOfSBATCHscript.sh 

#--------------------------------------
# 2. Convert downloaded sra files to fastq files 
#--------------------------------------
# Files will need to be converted from sra to fastq for downstream analysis. FASTQ format is a text-based format for storing both a biological sequence (usually nucleotide sequence) 
# and its corresponding quality scores. Both the sequence letter and quality score are each encoded with a single ASCII character for brevity.
#
# Example commmand: $fastq-dump sampleID.sra
# fastq-dump                                                     program part of the SRA toolkit that converts SRA files to fastq format
# sampleID.sra                                                   path and name of sample.sra that you would like to convert to fastq format
#
# will output sampleID.fastq 
# .fastq files were saved in fastq_files directory in Project /Project/fastq_files

# SBATCH script - 
sbatch NameOfSBATCHscript.sh 

#--------------------------------------
# 3. Create and view fastqc reports
#--------------------------------------
# Fastqc reads raw sequence data from high throughput sequencers and runs a set of quality checks to produce a report. 
# Best reports are those whose "per base sequence quality" are included in the green area of the graph & kmer content is good or average.
#
# Example command: $fastqc sampleID.fastq
# fastqc                                                        Babraham bioinformatics program that that checks for quality of reads 
# sampleID.fastq                                                path and name of sampleID in fastq format, may also be in fastq.gz format
#
# Reports were saved in fastq_files directory in Project /Project/fastq_files
# This command will create two outputs: an .html file & an .zip file. Will output sampleID_fastqc.html and sampleID_fastqc.zip files

# SBATCH script -
sbatch NameOfSBATCHscript.sh 

# Move fastqc reports to desktop to visualize them as you can't open html in a terminal.
# Open new terminal as this will not work if logged into a HPC (high performance computing) cluster
# Example command: $scp user@saguaro.a2c2.asu.edu:/Project/fastq_files/sampleID_raw_fastqc.html /Users/Desktop/
# scp                                                           secure copy  (linux command)                   
# /path/to/fastqc.html                                          path to where the files are located
# /path/where/you/want/fastqc.html                              path to where you would like to copy the files to 
#
# Reports were saved in Desktop in folder /raw_FASTQC
#--------------------------------------
# 4. Trim fastq files for quality and to remove adaptors 
#--------------------------------------
# In this case, fastq files need to be trimmed because the kmer count is poor, meaning that the RNA-seq has too many adapter sequences that need to be cut.
# To correct for this, trim the reads using Trimmomatic. Trimmomatic performs a variety of useful trimming tasks for illumina paired-end and single ended data.T
# The selection of trimming steps and their associated parameters are supplied on the command line.
# 
# For single-ended data, one input and one output file are specified, plus the filtering options. For paired-end data, two input files are specified (one file for each pair-end), and 4 output files, 2 for the 'paired' output where both reads survived the processing, and 2 for corresponding 'unpaired' output where a read survived, but the partner read did not.
# For this project, we have single-ended data.
#
# The current trimming steps are:
# ILLUMINACLIP: Cut adapter and other illumina-specific sequences from the read.
# SLIDINGWINDOW: Perform a sliding window trimming, cutting once the average quality within the window falls below a threshold.
# LEADING: Cut bases off the start of a read, if below a threshold quality
# TRAILING: Cut bases off the end of a read, if below a threshold quality
# CROP: Cut the read to a specified length
# HEADCROP: Cut the specified number of bases from the start of the read
# MINLEN: Drop the read if it is below a specified length
# TOPHRED33: Convert quality scores to Phred-33
# TOPHRED64: Convert quality scores to Phred-64
# It works with FASTQ (using phred + 33 or phred + 64 quality scores, depending on the Illumina pipeline used), either uncompressed or gzipp'ed FASTQ. Use of gzip format is determined based on the .gz extension.
# 
# The parameters selected were slidingwindow:4:30 leading10 trailing25 minlen50 phred33. They were chosen based on a better per sequence quality base and kmer content.
# A few fastq files were tested to determine the use of phred33 or phred 64. Phred33 was chosen due to the amount of bases ramining after trimming.
# Other parameters tested but not found to be ideal are sliding4:30 leading30 trailing40, sliding4:40 leading10 trailing25, sliding 4:40 leading30 trailing 40
# 
# Example Command: $java -jar trimmomatic-0.36.jar SE -phred33 input.fq.gz output.fq.gz ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:10 TRAILING:25 SLIDINGWINDOW:4:30 MINLEN:50
# java                                                          indicates that this is a java program and will require java in order to run
# -jar                                                          jar file to follow
# trimmomatic-0.36.jar                                          tool that will trim the raw fastq files
# SE                                                            SE is for singel end reads. If pair end then PE
# -phred33                                                      Using phred + 33 or phred + 64 quality scores, depending on the Illumina pipeline used, either uncompressed or gzipp'ed FASTQ 
# input.fq.gz                                                   sampeID in fastq format
# output.fq.gz                                                  sampleID output file. Use a descriptive name such as sampleID_minlen50_sliding430_leading30_trailing40.fq
# ILLUMINACLIP:TruSeq3-SE:2:30:10                               Remove Illumina adapters provided in the TruSeq3-PE.fa file (provided). Initially Trimmomatic will look for seed matches (16 bases) allowing maximally 2 mismatches. These seeds will be extended and clipped if in the case of paired end reads a score of 30 is reached (about 50 bases), or in the case of single ended reads a score of 10, (about 17 bases).
# LEADING:10                                                    Cut bases off the start of a read, if below a threshold quality of 10
# TRAILING:25                                                   Cut bases off the end of a read, if below a threshold quality of 25
# SLIDINGWINDOW:4:30                                            Scan the read with a 4-base wide sliding window, cutting when the average quality per base drops below 30
# MINLEN:50                                                     Drop the read if it is below a specified length of 50

# SBATCH script - 
sbatch NameOfSBATCHscript.sh 

#--------------------------------------
# 5. Create and view fastqc reports for the trimmed fastq files
#--------------------------------------
# Fastqc reads sequence data from high throughput sequencers and runs a set of quality checks to produce a report. 
# Best reports are those whose "per base sequence quality" are included in the green area of the graph & kmer content is good or average.
#
# Example command: $fastqc sampleID_minlen50_sliding430_leading30_trailing40.fastq
# fastqc                                                        Babraham bioinformatics program that that checks for quality of reads 
# sampleID.fastq                                                path and name of sampleID in fastq format, may also be in fastq.gz format
#
# Reports were saved in fastq_files directory in Project /Project/fastq_files
# This command will create two outputs: an .html file & an .zip file. Will output sampleID_fastqc.html and sampleID_fastqc.zip files

# SBATCH script -
sbatch NameOfSBATCHscript.sh 

# Move fastqc reports to desktop to visualize them as you can't open html in a terminal.
# Open new terminal as this will not work if logged into a HPC (high performance computing) cluster
# Example command: $scp user@saguaro.a2c2.asu.edu:/Project/fastq_files/sampleID_raw_fastqc.html /Users/Desktop/
# scp                                                           secure copy  (linux command)                   
# /path/to/fastqc.html                                          path to where the files are located
# /path/where/you/want/fastqc.html                              path to where you would like to copy the files to 
#
# Reports were saved in Desktop in folder /trimmed_FASTQC
#--------------------------------------
# 6. Obtain reference genome and gene annotation file 
#--------------------------------------
# Obtain reference genome and gene annotation file to be used for mapping reads. 
# Use the most relavent current version GRCh38.p7 from gencode. http://www.gencodegenes.org/releases/current.html
# Genome sequence (GRCh38.p7)	ALL	GRCh38.p7.genome.fa.gz
# Nucleotide sequence of the GRCh38.p7 genome assembly version on all regions, including reference chromosomes, scaffolds, assembly patches and haplotypes
# The sequence region names are the same as in the GTF/GFF3 files
# Example command $wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/GRCh38.p7.genome.fa.gz
# wget                                                          stands for "web get"
# ftp                                                           stands for File transfer program
# /GRCh38.p7.genome.fa.gz										                    human reference genome in fasta format

# Comprehensive gene annotation	ALL	gencode.v25.chr_patch_hapl_scaff.annotation.gtf.gz
# It contains the comprehensive gene annotation on the reference chromosomes, scaffolds, assembly patches and alternate loci (haplotypes)
# This is a superset of the main annotation file
# Example command $wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/gencode.v25.chr_patch_hapl_scaff.annotation.gtf.gz
# wget                                                          stands for "web get"
# ftp                                                           stands for File transfer program
# gencode.v25.chr_patch_hapl_scaff.annotation.gtf.gz			      human reference gene annotation file in .gtf format 

# SBATCH script -
sbatch NameOfSBATCHscript.sh 

#--------------------------------------
# 7. Generate reference genome index and dictionary 
#--------------------------------------
# A .dict dictionary of the contig names and sizes and a .fai fasta index file allow efficient random access to the reference bases for downstream analysis and mapping 
# You have to generate these files in order to be able to use a Fasta file as reference.
# In this step user supplied the reference genome sequences (FASTA files) and annotations(GTF file), from which STAR generate genome indexes that are utilized in the 2nd step. 
# The genome indexes are saved to disk and need only be generated once for each genome/annotation combination. 
# Note - files can not be be zipped (i.e file.gz) use $gunzip file.gz 
# Example command: $gunzip reference.gz
# gunzip														                            linux command to unzip zipped files
# reference.gz													                        path and name to zipped file
#
# Index reference genome using STAR
# Example command: $STAR --runMode genomeGenerate --genomeDir /GRCh38.p7/ --genomeFastaFiles GRCh38.p7.genome.fa.gz --sjdbGTFfile gencode.v25.chr_patch_hapl_scaff.annotation.gtf.gz --runThreadN 8
# STAR															                             path to calling STAR read aligner package
# --runMode 
# genomeGenerate												                         option directs STAR to run genome indices generation job.
# --genomeDir 
# /path/to/genomeDir
# --genomeFastaFiles 
# /path/to/genome/fasta 
# --sjdbGTFfile /path/to/annotations.gtf
# --sjdbOverhang ReadLength-1
# --runThreadN 													                         NumberOfThreads. Option defines the number of threads to be used for genome generation, it has to be set to the number of available cores on the server node.

# Create reference dictionrary using Picard tools CreateSequenceDictionary 
# Example command: $java -Xmx14g -jar picard.jar CreateSequenceDictionary R=GRCh38.p7.genome.fa O=GRCh38.p7.genome.fa.dict
# java 															                            java program and will need java installed in order to run
# -Xmx14g 														                          declares memory 
# -jar 															                            tells the program that the following tool is in .jar format
# picard.jar 													                          picard tools from the Broad Institute. A set of command line tools (in Java) for manipulating high-throughput sequencing (HTS) data and formats such as SAM/BAM/CRAM and VCF.
# CreateSequenceDictionary										                  creates a sequence dictionary for a reference sequence. The output file contains a header but no SAMRecords, and the header contains only sequence records.
# R=GRCh38.p7.genome.fa 										                    R stands for reference. Path and name of reference genome 
# O=GRCh38.p7.genome.fa.dict									                  O stand for output. Path and name of the output reference genome dictionary

# SBATCH script -
sbatch NameOfSBATCHscript.sh 

#--------------------------------------
# 8. STAR: map transcript reads to the reference genome, output as a .bam
#--------------------------------------
# STAR read aligner is a 2 pass process. 
# The user supplies the genome files generated in the pervious step (generate genome indexes), as well as the RNA-seq reads (sequences) in the form of FASTA or FASTQ files. 
# STAR maps the reads to the genome, and writes several output files, such as alignments (SAM/BAM), mapping summary statistics, splice junctions, unmapped reads, signal (wiggle) tracks etc. 
# Mapping is controlled by a variety of input parameters (options)

# STAR highly recommends using --sjdbGTFfile which specifies the path to the file with annotated transcripts in the standard GTF format. 
# Where STAR will extract splice junctions from this file and use them to greatly improve accuracy of the mapping. 
# While this is optional, and STAR can be run without annotations, using annotations is highly recommended whenever they are available.
# However this option should not be included for projects that include hybrids, as this might cause a bias towards the reference. 

# Compatibility with Cufflinks/Cuffdiff.
# For unstranded RNA-seq data, Cufflinks/Cuffdiff require spliced alignments with XS strand attribute, which STAR will generate with --outSAMstrandField intronMotif option. 
# As required, the XS strand attribute will be generated for all alignments that contain splice junctions. The spliced alignments that have undefined strand (i.e. containing only non-canonical unannotated junctions) will be suppressed.
# If you have stranded RNA-seq data, you do not need to use any specific STAR options. Instead, you need to run Cufflinks with the library option --library-type options. 
# For example, cufflinks ... --library-type fr-firststrand should be used for the standard dUTP protocol, including Illumina’s stranded Tru-Seq. This option has to be used only for Cufflinks runs and not for STAR runs.
# In addition, it is recommended to remove the non-canonical junctions for Cufflinks runs using --outFilterIntronMotifs RemoveNoncanonical.

# Unstranded vs. stranded RNAseq
# stranded vs. unstranded RNAseq data is where the strand specificity of origin for each transcript is defined. 
# Without strand information it is difficult and sometimes impossible to accurately quantify gene expression levels 
# for genes with overlapping genomic loci that are transcribed from opposite strands (Zhao et al. 2015)
# Stranded RNA-seq provides a more accurate estimate of transcript expression compared with non-stranded RNA-seq, 
# and is therefore the recommended RNA-seq approach for future mRNA-seq studies.
# If your data is unstranded then you will want to include the STAR options --outSAMstrandField intronMotif option.
# STAR will generate the XS strand attribute for all alignments that contain splice junctions. 
# The spliced alignments that have undefined strand (i.e. containing only non-canonical unannotated junctions) will be suppressed. 
# If you have stranded RNA-seq data, you do not need to use any specific STAR options. Instead, you need to run Cufflinks with the library option --library-type options (Dobin et al. 2012)

# First pass: maps fastq files to the reference genome and identifies splice junctions.
# Example comamnd: $STAR --genomeDir Project_genome --genomeLoad LoadAndKeep --readFilesIn sampleID.fastq --outSAMtype BAM Unsorted --outFileNamePrefix sampleID_pass1. --runThreadN 8
# STAR 									                                      path to calling STAR read aligner package
# --genomeDir 							                                  Specifies path to the directory (henceforth called ”genome directory” where the genome indices are stored. This directory has to be created (with mkdir) before STAR run and needs to have writing permissions. 
# Project_genome												                      path and directory to reference genome
# --genomeLoad													                      mode of shared memory usage for the genome files
# LoadAndKeep													                        load genome into shared and keep it in memory after run
# --sjdbGTFfile													                      path to the GTF file with annotations
# --readFilesIn 						                                  if pair end reads, include path to both reads with a " " space inbetween _1 _2, /geuvadis_fastq/sampleID_1.fastq /home/kcolney/map_geuvadis/geuvadis_fastq/sampleID_2.fastq 
# sampleID.fastq												                      name and path to sample in fastq format. If paired end samples include both pairs and separate with a space i.e (sample_1.fastq sample_2.fastq)
# --outSAMtype 							                                  indicate which output format, BAM unsorted
# BAM Unsorted													                      output unsorted Aligned.out.bam file. The paired ends of an alignment are always adjacent, and multiple alignments of a read are adjacent as well. This ”unsorted” file can be directly used with downstream software such as HTseq, without the need of name sorting. The order of the reads will match that of the input FASTQ(A) files only if one thread is used
# --outSAMstrandField 											                  * This option depends on the data * For unstranded RNA-seq data, Cufflinks/Cuffdiff require spliced alignments with XS strand attribute, which STAR will generate with --outSAMstrandField intronMotif option
# intronMotif													                        * This option depends on the data * As required, the XS strand attribute will be generated for all alignments that contain splice junctions. The spliced alignments that have undefined strand (i.e. containing only non-canonical unannotated junctions) will be suppressed.
# --outFilterIntronMotifs 										                * This option depends on the data * it is recommended to remove the non-canonical junctions for Cufflinks runs using
# RemoveNoncanonical											                    * This option depends on the data * remove the non-canonical junctions
# --outFileNamePrefix 					                              define the sample id prefix, sampleID_pass1. (bam, will be added by the STAR program)
# --runThreadN 							                                  for computing purpose allocate the number of threads, 8 

# SBATCH script - 
sbatch NameOfSBATCHscript.sh 

# For the most sensitive novel junction discovery, it is recommend to run STAR in the 2-pass mode. It does not increase the number of detected novel junctions, but allows to detect more splices reads mapping to novel junctions. 
# The basic idea is to run 1st pass of STAR mapping with the usual parameters, then collect the junctions detected in the first pass, and use them as ”annotated” junctions for the 2nd pass mapping.
# Second pass: maps fastq files between each other and identifies splice junctions.
# STAR will utilize annotations formatted as a list of splice junctions coordinates in a text file: --sjdbFileChrStartEnd /path/to/sjdbFile.txt. 
# This file should contains 4 columns separated by tabs:
# Chr \tab Start \tab End \tab Strand=+/-/.
# Example command: $STAR --genomeDir Project_genome --readFilesIn sampleID_1.fastq --outSAMtype BAM Unsorted --outFileNamePrefix sampleID_pass2. --sjdbFileChrStartEnd /sampleID1_pass1.SJ.out.tab /sampleID2_pass1.SJ.out.tab /sampleID3_pass1.SJ.out.tab --runThreadN 14
# STAR 									                                  STAR read aligner package
# --genomeDir 							                              define where the genome is location, /star_genome 
# Project_genome												                  path and directory to reference genome
# --genomeLoad													                  mode of shared memory usage for the genome files
# LoadAndKeep													                    load genome into shared and keep it in memory after run
# --sjdbGTFfile													                  path to the GTF file with annotations
# --readFilesIn 						                              if pair end reads, include path to both reads with a " " space inbetween _1 _2, /geuvadis_fastq/sampleID_1.fastq /home/kcolney/map_geuvadis/geuvadis_fastq/sampleID_2.fastq 
# sampleID_1.fastq												                name and path to sample in fastq format. If paired end samples include both pairs and separate with a space i.e (sample_1.fastq sample_2.fastq)
# --outSAMtype 							                              indicate which output format, BAM unsorted
# BAM Unsorted													                  output unsorted Aligned.out.bam file. The paired ends of an alignment are always adjacent, and multiple alignments of a read are adjacent as well. This ”unsorted” file can be directly used with downstream software such as HTseq, without the need of name sorting. The order of the reads will match that of the input FASTQ(A) files only if one thread is used
# --sjdbFileChrStartEnd					                          path to the pass_1.SJ.out.tab files made in the first pass
# sampleID1_pass1.SJ.out.tab									            4 columns separated by tabs: Chr \tab Start \tab End \tab Strand=+/-/. Here Start and End are first and last bases of the introns (1-based chromosome coordinates). This file can be used in addition to the --sjdbGTFfile, in which case STAR will extract junctions from both files.
# sampleID2_pass1.SJ.out.tab									            List all the samples from the first pass or all the samples in a group (i.e population, cases and controls, hybrids, males and females)
# --outFileNamePrefix 					                          define the sample id prefix, sampleID_pass1. (bam, will be added by the STAR program)
# --runThreadN 							                              for computing purpose allocate the number of threads, 14 

# SBATCH script - 
sbatch NameOfSBATCHscript.sh 

#--------------------------------------
# 9. Check initial quality stats on .bam files 
#--------------------------------------
# Will print basic statistics from input BAM file(s)
# Total reads:       
# Mapped reads:      
# Forward strand:    
# Reverse strand:    
# Failed QC:         
# Duplicates:        
# Paired-end reads:  
# 'Proper-pairs':    
# Both pairs mapped: 
# Read 1:            
# Read 2:            
# Singletons:   
#     
# Example command: $bamtools stats -in sampleID_pass2.Aligned.out.bam > sampleID_pass2.txt
# bamtools											                      package 
# stats												                        command to get general alignment statistics
# -in												                          indicates input file
# sampleID.bam										                    path and name to bam file 
# >												                            directs output     
# sampleID_pass2.txt                                  indicated output file name

# SBATCH script - 
sbatch NameOfSBATCHscript.sh 

#--------------------------------------
# 10. Sort bam files
#--------------------------------------
# An appropriate @HD-SO sort order header tag will be added or an existing one updated if necessary.
# Example command: $bamtools sort -in sampleID.bam -out sampleID.sorted.bam
# bamtools								                          package 
# sort								                              command to add header tags
# -in									                              indicates input file
# sampleID.bam						                          path and name to bam file 
# -out 									                            indicated output file
# sampleID.sorted.bam					                      the sorted BAM file 

# SBATCH script - 
sbatch NameOfSBATCHscript.sh 

#--------------------------------------
# 11. Generate stats on sorted bam files
#--------------------------------------
# Will print basic statistics from input BAM file(s)   
# Example command: $bamtools stats -in sampleID_pass2.Aligned.sorted.bam > sampleID_pass2.Aligned.sorted.txt
# bamtools											                    package 
# stats												                      command to get general alignment statistics
# -in												                        indicates input file
# sampleID.bam										                  path and name to bam file 
# >												                          directs output     
# sampleID_pass2.txt                                indicated output file name

# SBATCH script - 
sbatch NameOfSBATCHscript.sh 

# Compare sorted.bam stats to the original .bam stats, there should be no differences between them
# We do this to step (bamtools stats) every time we do anything to our bam files as a quality control check
#--------------------------------------
# 12. Mark duplicates & Remove duplicates
#--------------------------------------
# Mark duplicates: "Flags" where the duplicate reads are
# Example command: $java -Xmx8g -jar picard.jar MarkDuplicates INPUT=sampleID.sorted.bam OUTPUT=sampleID.sorted.markdup.bam METRICS_FILE=sampleID.markdup.picardMetrics.txt REMOVE_DUPLICATES=false ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT
# java												                      program called 
# -Xmx8g 											                      declares memory 
# picard.jar 										                    path to picard jar file 
# MarkDuplicates 									                  command to create sequence dictionary 
# INPUT=											                      path to input file, sorted bam file per sample	
# OUTPUT= 											                    path and name or output file .markdup to indicate this file will contain duplicates that have been marked
# METRICS_FILE=										                  file to write duplication metrics to save as sampleID.markdup.picardMetrics.txt
# REMOVE_DUPLICATES=false 							            If true do not write duplicates to the output file instead of writing them with appropriate flags set. Default value: false. This option can be set to 'null' to clear the default value. Possible values: {true, false}
# ASSUME_SORTED=true 								                BAM files are sorted because we sorted them in step 6
# VALIDATION_STRINGENCY=LENIENT						          setting stringency to SILENT can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded.

# SBATCH script - 
sbatch NameOfSBATCHscript.sh 

#--------------------------------------
# 13. Generate statistics for bam files with duplicates marked 
#--------------------------------------
# For each sample get the read stats for the mark duplicates BAM files 
# Example command: $bamtools stats -in sampleID_pass2.Aligned.sorted.dupsMark.bam > sampleID_pass2.Aligned.sorted.dupsMark.txt

# SBATCH script - 
sbatch NameOfSBATCHscript.sh 

# Compare stat results for each sample to the original bam file (sanity check: is there the same number of reads as the original BAM file?)
# If there is more than 15% of the reads being marked as duplicates may need to consider removing that sample
#--------------------------------------
# 14. Add or replace read groups
#--------------------------------------
# For each sample, add a read group to the mark duplicate BAM files (a read group is a "tag" such as a sample ID)
# Example command: $java -Xmx8g -jar picard.jar AddOrReplaceReadGroups INPUT=sampleID.sorted.rmdup.bam OUTPUT=sampleID.sorted.rmdup.addReadGr.bam RGLB=sampleID RGPL=machineUsed RGPU=laneUsed RGSM=sampleName RGCN=location RGDS=species VALIDATION_STRINGENCY=LENIENT
# java												                      program called 
# -Xmx8g 											                      declares memory 
# picard.jar 										                    path to picard jar file 
# AddOrReplaceReadGroups							              Replaces all read groups in the INPUT file with a single new read group and assigns all reads to this read group in the OUTPUT BAM
# INPUT=											                      path to input file, sorted bam file per sample	
# OUTPUT= 											                    path and name or output file .markdup to indicate this file will contain duplicates that have been marked
# RGLB=												                      Read Group Library Required (sampleID)
# RGPL=												                      Read Group platform (e.g. illumina, solid) Required 
# RGPU=	 											                      Read Group platform unit (eg. run barcode) Required	(laneUsed)
# RGSM=	 											                      Read Group sample name Required (sampleName or sampleID)
# RGCN= 											                      Read Group sequencing center name Default value: null (i.e. ASU)
# RGDS=												                      Read Group description Default value: null (speciesName)
# VALIDATION_STRINGENCY=LENIENT						          setting stringency to SILENT can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded.

# SBATCH script - 
sbatch NameOfSBATCHscript.sh 

#--------------------------------------
# 15. Generate statistics on read group BAM sample files
#--------------------------------------
# For each sample get the read stats for the remove duplicates and add read groups BAM files.
# Statistics on the BAM files should be the same as before the previous step when read groups were modified.
# Compare stat results for each sample to the markdup.bam file (sanity check: is there the same number of reads as the original BAM file?)
# Example command: bamtools stats -in sampleID.sorted.markdup.addReadGr.bam

# SBATCH script - 
sbatch NameOfSBATCHscript.sh 

#--------------------------------------
# 16. Index BAM files 
#--------------------------------------
# For each sample index the processed BAM files that are sorted, have marked duplicates, and have read groups added. These will be used to identify callable loci. 
# Indexing is used to "sort" by chromosome and region 
# Output will be sampleID.sorted.markdup.addReadGr.bam.bai
# Example command: $bamtools index -in sampleID.sorted.markdup.addReadGr.bam
# bamtools											                    package 
# index												                      Generates index for BAM file
# -in												                        indicates input file
# sampleID.sorted.markdup.addReadGrbam				      path and name to bam file that has been sorted, duplicate reads were removed, and the read groups were added 

# SBATCH script - 
sbatch NameOfSBATCHscript.sh 

#--------------------------------------
# 17. SplitNCigarReads - splits reads into exon segments 
#--------------------------------------
# SplitNCigarReads splits reads into exon segments (getting rid of Ns but maintaining grouping information) and hard-clip any sequences overhanging into the intronic regions. 
# Overhang regions that have mismatch greater than the threshold are clipped to reduce the number of the called false variants. 
# At SplitNCigarReads, mapping qualities are also reassigned because STAR assigns good alignments a MAPQ of 255 which is meaningless to GATK. 
# To correct for this, GATK’s ReassignOneMappingQuality read filter, is used to reassign all good alignments to the default value of 60. 
# Example command: $java -Xmx24g -Djava.io.tmpdir=~/temp/ -jar /GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar -T SplitNCigarReads -R /reference/referenceGenome.fa -I sampleID.sorted.markdup.addReadGr.bam -o sampleID.sorted.markdup.addReadGr.split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS
# java												                                        program called 
# -Xmx8g 											                                        declares memory 
# -Djava.io.tmpdir=~/temp/ 							                              Temporary directory for extra storage  
# -jar 												                                        declares jar file type 
# /GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar 		                      GATK's GenomeAnalysisTK.jar tool 
# -T 												                                          declares tool 
# SplitNCigarReads 									                                  splits reads into exon segments (getting rid of Ns but maintaining grouping information) and hard-clip any sequences overhanging into the intronic regions
# -R 												                                          reference 
# /reference/referenceGenome.fa 					                            path and genome reference file 
# -I 												                                          input 
# sampleID_pass2.sorted.rg.dupsMarked.bam 			                      input bam file, for each sample, the pass2. sorted. read groups add and duplicates have been marked. 	
# -o												                                          output	
# sampleID_pass2.sorted.rg.dupsMarked.split.bam 	                    clearly identify the output file name. include the sample ID and .bam 	
# -rf 												                                        declares tool 	
# ReassignOneMappingQuality 						                              ReassignOneMappingQuality read filter to reassign all good alignments to the default value of 60
# -RMQF 255 										                                      STAR assigns good alignments a MAPQ of 255 
# -RMQT 60 											                                      reassign all good alignments to the default value of 60 	
# -U ALLOW_N_CIGAR_READS							                                reads with N cigars should be allowed	

# SBATCH script - 
sbatch NameOfSBATCHscript.sh 

#---------------------------------------
# 18. Differential expression using Cuffdiff
#---------------------------------------
# Identify differentially expressed genes and transcripts using cuffdiff. 
# Cuffdiff is used to find significant changes in transcript expression, splicing, and promoter use. 
# Cuffdiff uses a gene annotation file downloaded with a sbatch script from UCSC hg19
# Example command: $cuffdiff -use-sample-sheet -o diff_out -b reference.fa -p 8 --library-type fr-firststrand -L set1,set2 -u refence.gtf sampleSet.txt
# cuffdiff 											                                    program in cufflinks package that identifies if genes are differentially expressed 
# -use-sample-sheet 								                                tells the program to use a text file containing the path to the samples and the group_id 
# -o 												                                        output
# diff_out 											                                    path and name of output directory where all of the output files from the cuffdiff program will be placed
# -b 												                                        Providing Cufflinks with the multifasta file your reads were mapped to via this option instructs it to run the bias detection and correction algorithm which can significantly improve accuracy of transcript abundance estimates
# reference.fa 										                                  path and name of reference.fasta which the files were mapped to
# -p 												                                        Use this many threads to align reads. The default is 1
# 8 												                                        declared 8 threads 
# --library-type												                            * This option depends on the data * unstranded RNA-seq data, Cufflinks/Cuffdiff require spliced alignments with XS strand attribute
# fr-firststrand												                            * This option depends on the data * should be used for the standard dUTP protocol, including Illumina’s stranded Tru-Seq.
# -L 												                                        Specify a label for each sample, which will be included in various output files produced by Cuffdiff.
# set1,set2 										                                    Must have at least 2 labels. (i.e male,female)	
# -u 												                                        Tells Cufflinks to do an initial estimation procedure to more accurately weight reads mapping to multiple locations in the genome
# refence.gtf 										                                  path and name to the gene annotation file. 
# sampleSet.txt										                                  text file containing the list of samples

# SBATCH script - 
sbatch NameOfSBATCHscript.sh 

#---------------------------------------
# 19. Variant calling using GATK's HaplotypeCaller
#---------------------------------------
# Call germline SNPs and indels via local re-assembly of haplotypes
# Variant calling is  accomplished using HaplotypeCaller tool. The HaplotypeCaller tool takes into account the information about intron-exon split regions that is embedded in the BAM file by SplitNCigarReads (pervious step). 
# HaplotypeCaller is capable of calling SNPs and indels simultaneously via local de-novo assembly of haplotypes in an active region.
# HaplotypeCaller runs per-sample to generate an intermediate genomic gVCF (gVCF), which can then be used for joint genotyping of multiple samples (next step). 
# HaplotypeCaller Steps:
#		1. Defines active regions, ased on evidence for variation determines which regions in the genome it needs to perform on. 
#		2. Determines haplotypes by assembly of the active region. For each ActiveRegion, the program builds a De Bruijn-like graph to reassemble the ActiveRegion, and identifies what are the possible haplotypes present in the data. 
#		   The program then realigns each haplotype against the reference haplotype using the Smith-Waterman algorithm in order to identify potentially variant sites.
#		3. Determines the likelihoods of the haplotypes given the read data. Performs a pairwise alignment of each read against each haplotype using the PairHMM algorithm. 
#		   Produces a matrix of likelihoods of haplotypes given the read data. These likelihoods are then marginalized to obtain the likelihoods of alleles for each potentially variant site given the read data.
#		4. Assigns samples genotypes. For each potentially variant site, the program applies Bayes' rule, using the likelihoods of alleles given the read data to calculate the likelihoods of each genotype per sample given the read data observed for that sample. 
#		   The most likely genotype is then assigned to the sample.
# Example command: $java -Xmx24g -Djava.io.tmpdir=~/temp/ -jar /GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar -T HaplotypeCaller -R /reference/referenceGenome.fa -I sampleID.sorted.markdup.addReadGr.split.bam -stand_call_conf 30.0 -stand_emit_conf 30.0  -o sampleID.gvcf
# java												                                      program called 
# -Xmx8g 											                                      declares memory 
# -Djava.io.tmpdir=~/temp/ 							                            Temporary directory for extra storage  
# -jar 												                                      declares jar file type 
# /GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar 		                    GATK's GenomeAnalysisTK.jar tool 
# -T 												                                        declares tool 
# HaplotypeCaller 									                                Calls SNPs and indels simultaneously via local de-novo assembly of haplotypes in an active region
# -R 												                                        reference
# /reference/referenceGenome.fa 					                          path and genome reference file
# -I 												                                        input
# sampleID_pass2.sorted.rg.dupsMarked.split.bam  	                  input bam file, for each sample, the pass2. sorted. read groups added, duplicates have been marked, and SplitNCigarReads has been preformed
# -stand_call_conf 30.0 							                              The minimum phred-scaled confidence threshold at which variants should be called
# -stand_emit_conf 30.0 							                              The minimum phred-scaled confidence threshold at which variants should be emitted (and filtered with LowQual if less than the calling threshold)
# --emitRefConfidence GVCF 							                            Mode for emitting reference confidence scores
# -o                                                                indicates output 
# sampleID.gvcf                                                     path and name of sample and followed by .gvcf

# SBATCH script - 
sbatch NameOfSBATCHscript.sh 

#----------------------------------------
# 20. Merge gvcf files to a group VCF file 
#----------------------------------------
# GVCF stands for Genomic VCF. A GVCF is a kind of VCF, so the basic format specification is the same as for a regular VCF, but a Genomic VCF contains extra information.
# Merge the gvcf sample files into one raw vcf file. 
# Input the --variant sampleID.gvcf in the order in which you want them in your raw vcf output file
# Example command: $java -Xmx16g -jar GenomeAnalysisTK.jar -T GenotypeGVCFs -R referenceGenome.fna --variant sampleID.rmdup.gvcf --variant sampleID.rmdup.gvcf  -o Project.raw.vcf
# java												                                    program called 
# -Xmx16g 											                                  declares memory
# GenomeAnalysisTK.jar 								                            path to GATK jar file 
# -T												                                      declares that a tool command will follow 
# GenotypeGVCFs										                                GenotypeGVCFs merges gVCF records that were produced in the previous step. This tool performs the multi-sample joint aggregation step and merges the records together. 
#													                                        At each position of the input gVCFs, this tool will combine all spanning records, produce correct genotype likelihoods, re-genotype the newly merged record, and then re-annotate it.	
# -R 												                                      reference genome will follow this tag
# referenceGenome									                                path and name of reference genome 
# --variant											                                  One or more input gVCF files. Must be declared before each sample input and then followed by the sample gvcf file 
# sampleID.rmdup.gvcf								                              sample gvcf file made in pervious step
# -o												                                      declares output to follow 
# Project.raw.vcf								                                  Name of output. Include the .vcf and include "raw" because this is an unfiltered vcf file. 

# SBATCH script - 
sbatch NameOfSBATCHscript.sh 

#----------------------------------------
# 21. Remove indels and only keep biallelic sites from VCF - BCF file 
#----------------------------------------
# An indel is the insertion or deletion of bases in the DNA of an organism. 
# It has slightly different definitions between its use in evolutionary studies and its use in germ-line and somatic mutation studies.
# Indels will need to be filtered out of the merged.vcf because when we later use GATK's ASEReadCounter for obtaining counts per allele for allele-specific expression analysis
# will not be able to properly handle indels and will call them as being homozygous at that site. Therefore, they will need to be filtered out using BCFtools. 
#
# -v, --types snps|indels|mnps|other
# comma-separated list of variant types to select. Site is selected if any of the ALT alleles is of the type requested. 
# 
# Example command: $bcftools view Project.raw.vcf -m2 -M2 -v snps > Project.ONLYbiallelicSNPSites.vcf
# bcftools											                              program called 
# view 												                                tells bcftools to open and view the the VCF that needs to be filtered
# Project.raw.vcf						                                  name of the raw VCF 
# -m2												                                  prints at most 2 alleles in the REF column (reference) 
# -M2												                                  prints at most 2 alleles in the ALT column (alternate)
# -v 												                                  selects for "types". Types are determined by comparing the REF and ALT alleles in the VCF record not INFO tags like INFO/INDEL or INFO/VT
# snps 												                                type snps selects for only SNP sites and will exclude all indels 
# >													                                  directs to output 
# Project.ONLYbiallelicSNPSites.vcf		                    	  name and path of output VCF that now only contains biallelic sites and no indels 

# SBATCH script - 
sbatch NameOfSBATCHscript.sh 

#----------------------------------------
# 22. Create 1 VCF per sample - SelectVariants 
#----------------------------------------
# Create 1 VCF per sample to help reduce how long it will take for GATK's ASEReadCounter to run 
# To do this, use GATK's SelectVariants. 
# Example command: $java -Xmx16g -jar GenomeAnalysisTK.jar -T SelectVariants -R gencode.GRCh38.fasta -V Project.ONLYbiallelicSNPSitess.vcf -o sampleID.ONLYbiallelicSNPSites.vcf -sn sampleID
# java												                                  program called 
# -Xmx16g 											                                declares memory
# -jar												                                  jar file 
# GenomeAnalysisTK.jar 								                          path to GATK jar file 
# -T 												                                    declares that a tool command will follow 
# SelectVariants 									                              Select a subset of variants from a larger callset
# -R												                                    reference genome to follow
# gencode.GRCh38.fasta 						                              path and name of reference genome 
# -V 												                                    A variant call set from which to select a subset
# Project.ONLYbiallelicSNPSitess.vcf 			                      path and name of the variant call set 
# -o 												                                    output to follow
# sampleID.ONLYbiallelicSNPSitess.vcf 		                      name of the new VCF file containing the selected subset of variants
# -sn 												                                  sample name select
# sampleID										                                  name of individual to be selected 

# SBATCH script - 
sbatch NameOfSBATCHscript.sh 

#----------------------------------------
# 23. Filter for only heterozygous sites - SnpSift
#----------------------------------------
# Using SnpSift filter you can filter VCF files using arbitrary expressions.
# We want to filter to only keep heterozygous sites because we will then run GATK's ASEReadCounter to count per allele for allele-specific expression analysis. 
# Example command: $cat sampleID.ONLYbiallelicSNPSitesss.vcf | java -jar SnpSift.jar filter "((countHet() > 0))" > sampleID.ONLYbiallelicSNP_HET_Sites.vcf
# cat                                                           cat stands for "catenate." It reads data from files, and outputs their contents.
# sampleID.ONLYbiallelicSNPSitesss.vcf                          path and name of sample.vcf
# |                                                             pipe the information generated in front of the "|" to be feed into the command following the pipe ("|")
# java                                                          java program called
# -jar                                                          file type is jar
# SnpSift.jar                                                   SnpEff program to filter a vcf
# filter                                                        filter tool
# "((countHet() > 0))"                                          keep only sites that are heterozygous
# >                                                             direct this to an output 
# sampleID.ONLYbiallelicSNP_HET_Sites.vcf                       name and path to output vcf that now will only contain heterozygous biallelic SNPs

# SBATCH script - 
sbatch NameOfSBATCHscript.sh 

#----------------------------------------
# 24. Allele Specific Expression 
#----------------------------------------
# GATK ASEReadCounter calculates allele counts at a set of positions after applying filters that are tuned for enabling allele-specific expression (ASE) analysis
# The filters operate on mapping quality, base quality, depth of coverage, overlapping paired reads and deletions overlapping the position. 
# This tool will only process biallelic SNP sites. If your callset contains multiallelic sites, they will be ignored
# ASEReadCounter only considers heterozygous sites and the minimum mapping quality default is 0 and the minimum base quality default is 0 
# Like most GATK tools, this tools filters out duplicate reads by default. However, some ASE methods recommend including duplicate reads in the analysis, so the DuplicateRead filter can be disabled using the "-drf DuplicateRead" flag in the command-line.
#
# Inputs:
# 	BAM files to be analyzed for ASE
# 	A VCF file with specific sites to process 
# 	Reference genome in fasta format
# Output:
# 	Table of allele counts at the given sites. Tab-delimited that is readable by R and compatible with Mamba (downstream tool developed for allele-specific expression analysis)
# Example command: $java -jar GenomeAnalysisTK.jar -R reference.fasta -T ASEReadCounter -o file_name.csv -I input.bam -sites sites.vcf -U ALLOW_N_CIGAR_READS -minDepth 10 --minMappingQuality 10 --minBaseQuality 2 -drf DuplicateRead
# java																																		program called
# -jar																																		declares jar file to follow
# GenomeAnalysisTK.jar 								            												path to GATK jar file
# -R 												            																	declares reference genome
# reference.fasta 									            													path to reference genome in fasta format
# -T 												            																	tool
# ASEReadCounter 									            														this tools filters out duplicate reads by default. However, some ASE methods recommend including duplicate reads in the analysis, so the DuplicateRead filter can be disabled using the "-drf DuplicateRead" flag in the command-line.
# -o 												            																	output 
# file_name.csv 									          															name of output file for each sample in csv format (must include .csv)
# -I 												            																	input
# input.bam 										          			 													input BAM file to be used 
# -sites 											            																declares vcf file to follow
# sites.vcf 										          																vcf file that was generated in pervious step 
# -U 												            																	Filtering options:
# ALLOW_N_CIGAR_READS 																										allow N cigar reads
# -minDepth 10                                                 						minimum depth
# --minMappingQuality 10																									minimum read mapping quality
# --minBaseQuality 2 																											minimum base quality
# -drf DuplicateRead																											like most GATK tools, this tools filters out duplicate reads by default. DuplicateRead filter can be disabled using the "-drf DuplicateRead" flag in the command-line.

# SBATCH script - 
sbatch NameOfSBATCHscript.sh 


