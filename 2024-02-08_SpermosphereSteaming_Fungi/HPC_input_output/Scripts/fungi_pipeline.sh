#!/bin/bash
#

##################################
#     Bioinformatic pipeline     #
#  Fungi ITS amplicon sequences  #
#                                #
#     ----------------------     #
#  creating a samples.txt file   #
#                                #
##################################

# It sets the variable PATH_SAMPLES_LIST to the path "../2023-11-09_QUO1007226_Bacteria_NCST2023_SoySeedEpiphytes".
# It uses the ls command to list the contents of the directory specified by PATH_SAMPLES_LIST.
# It pipes the output of ls to a series of sed (stream editor) commands:
# The first sed command "s/_R1_001.fastq.gz//" removes the "_R1_001.fastq.gz" substring from each line of the output.
# The second sed command "s/_R2_001.fastq.gz//" removes the "_R2_001.fastq.gz" substring from each line of the output.
# The resulting output, after the modifications by the sed commands, is a list of file names without the specified substrings.
# The uniq command is then used to filter out duplicate lines in the modified output.
# Finally, the unique list is redirected to a file named "samples.txt" using the > operator.
# In summary, this script generates a list of unique sample names by removing specific substrings from the file names in the specified directory and saves the unique list to a file named "samples.txt".

# Path to directory containing the samples you want in the samples.txt 
PATH_SAMPLES_LIST="steaming_fungi/"

echo "Path to samples is: $PATH_SAMPLES_LIST"
date

ls $PATH_SAMPLES_LIST | sed "s/_R1_001.fastq.gz//" | sed "s/_R2_001.fastq.gz//" | uniq > samples.txt

SAMPLES_LIST="samples.txt" 

##################################
#     Bioinformatic pipeline     #
#     ITS amplicon sequences     #
#                                #
#     ----------------------     #
#        stripping primers       #
#                                #
##################################


#It creates a directory named "trimmed" using the mkdir command in the current working directory (./ refers to the current directory).
#It sets a variable PATH_SAMPLES_LIST to the path of a file named "samples.txt".
#It loads the Anaconda module (specifically the Python 3 version from February 2020) using the module load command.
#It uses a for loop to iterate over each sample listed in the "samples.txt" file.
#Within the loop, it echoes the current sample being processed.
#It calls the cutadapt tool to perform adapter trimming on the forward fastq files for each sample:

mkdir ./trimmed

#  load the module
module load anaconda/3-2020.02

for sample in $(cat $SAMPLES_LIST)
do

    echo "On sample: $sample"

    cutadapt -g CTTGGTCATTTAGAGGAAGTAA -e 0.01 --discard-untrimmed --match-read-wildcards $PATH_SAMPLES_LIST/${sample}*.fastq.gz > trimmed/${sample}_trimmed.fastq

done

# Use this to generate a file to grab the number of passing reads
#grep -e "Reads written" -e "On sample:" 03strippingprimers.o1413 > percent_primers_passing.txt


##################################
#     Bioinformatic pipeline     #
#     ITS amplicon sequences     #
#                                #
#     ----------------------     #
#        stats                   #
#                                #
##################################

#It loads the "vsearch" module using the module load command.
#It creates a directory named "stats" using the mkdir command in the current working directory (./ refers to the current directory).
#It concatenates all the fastq files in the "trimmed" directory into a single file named "trimmed.fastq" using the cat command.
#It uses the "vsearch" tool to generate statistics for the concatenated fastq file:
#-fastq_stats trimmed/trimmed.fastq: Specifies the input fastq file for which statistics are generated.
#-log stats/stats_results.txt: Specifies the output file where the statistics results will be logged.
#The purpose of this script is to concatenate multiple trimmed fastq files into a single file and then use "vsearch" to generate statistics for the concatenated file. The statistics, such as the number of reads, mean quality scores, etc., will be written to the "stats_results.txt" file in the "stats" directory.

module load vsearch
mkdir ./stats

cat trimmed/*.fastq > trimmed/trimmed.fastq

vsearch -fastq_stats trimmed/trimmed.fastq -log stats/stats_results.txt


##################################
#     Bioinformatic pipeline     #
#     ITS amplicon sequences     #
#                                #
#     ----------------------     #
#        filtering               #
#                                #
##################################

# This script is designed to filter the "trimmed.fastq" file in the "trimmed" directory using vsearch, applying various filtering criteria. The filtered results are then saved in both FASTA and FASTQ formats in the "filtered" directory.

# This command loads the "vsearch" module, making the vsearch tool available for use in the script.
module load vsearch

#This command creates a directory named "filtered" in the current working directory (./ refers to the current directory).
mkdir ./filtered
mkdir ./trimmed

# Decisions on these options should be decided upon by the 04_stats script

#-fastq_filter trimmed/trimmed.fastq: Specifies the input fastq file to be filtered.
#-fastq_maxee 1: Filters out reads with an expected number of errors greater than 1.
#-fastq_trunclen 275: Truncates reads longer than 275 bases.
#-fastq_maxns 0: Filters out reads with more than 0 Ns (ambiguous bases).
#-fastaout filtered/filtered.fasta: Specifies the output file in FASTA format.
#-fastqout filtered/filtered.fastq: Specifies the output file in FASTQ format.

vsearch -fastq_filter trimmed/trimmed.fastq -fastq_maxee 1 -fastq_trunclen 264 -fastq_maxns 0 -fastaout filtered/filtered.fasta -fastqout filtered/filtered.fastq
vsearch -fastq_filter filtered/filtered.fastq -fastq_stripleft 44 -fastaout trimmed/trimmed_R1.fasta -fastqout trimmed/trimmed_R1.fastq


############################################
#     Bioinformatic pipeline                #
#     ITS amplicon sequences               #
#DEREPLICATION, CLUSTERING, CHIMERA REMOVAL#
############################################

# This script processes filtered sequences, dereplicates them, performs denoising, and clusters OTUs based on traditional 97% identity. 

source /apps/profiles/modules_asax.sh.dyn
module load usearch/11.0.667
module load vsearch

mkdir ./clustered

# --derep_fulllength: Performs dereplication of sequences.
# --output filtered/uniques.fasta: Specifies the output file for dereplicated sequences.
# -sizeout: Outputs sequence abundances in headers.
vsearch --derep_fulllength trimmed/trimmed_R1.fasta --output trimmed/uniques.fasta -sizeout


# de-noising (error correction), output is zero-radius OTUs

#-unoise3: Performs denoising using unoise3 algorithm.
#-tabbedout clustered/unoise_zotus_R1.txt: Outputs results in a tab-separated format.
#-zotus clustered/zotus_R1.fasta: Outputs zero-radius OTUs.

#usearch -unoise3 filtered/uniques.fasta -tabbedout clustered/unoise_zotus_R1.txt -zotus clustered/zotus_R1.fasta

# clusters OTUs based on traditional 97% identity 

#-cluster_otus: Performs OTU clustering.
#-minsize 2: Sets the minimum size of a cluster to 2.
#-otus clustered/otus.fasta: Specifies the output file for clustered OTUs.
#-uparseout filtered/uparse_otus.txt: Outputs UPARSE-specific information.
#-relabel BOTU_: Relabels OTUs with the prefix "BOTU_".
#--threads 20: Specifies the number of threads to use.

usearch -cluster_otus trimmed/uniques.fasta -minsize 2 -otus clustered/otus.fasta -uparseout trimmed/uparse_otus.txt -relabel FOTU_ --threads 20

# useful links
#http://www.drive5.com/usearch/manual/unoise_pipeline.html
#http://www.drive5.com/usearch/manual/faq_uparse_or_unoise.html
#http://www.drive5.com/usearch/manual/cmd_otutab.html
#http://www.drive5.com/usearch/manual/upp_labels_sample.html
#http://drive5.com/usearch/manual/bugs.html
#http://drive5.com/usearch/manual/support.html


##################################
#     Bioinformatic pipeline     #
#     		ITS Fungi            #
#     ----------------------     #
# taxonomy assignment - SINTAX   #
#                                #
#      		Zachary Noel         #
##################################

module load vsearch 
mkdir ./taxonomy

# Assign taxonomy using SINTAX algorithm
vsearch -sintax clustered/otus.fasta -db ~/noel_shared/db_fungi/sh_general_release_dynamic_s_all_25.07.2023_mockseqadded.fasta -tabbedout taxonomy/ITS_taxonomy.txt -strand both -sintax_cutoff 0.8

# Print the first and fourth columns of the sintax output, replace the commas with tabs, get rid of the taxonomic prefixes, save it in a text file
awk '{print $1, $4}' taxonomy/ITS_taxonomy.txt | tr ' ' ',' |  sed 's/d://; s/p://; s/c://; s/o://; s/f://; s/g://; s/s://' > taxonomy/ITS_taxonomy2.txt


############################################
#     Bioinformatic pipeline               #
#     	   NBC Taxonomy           		   #
#     									   #
############################################ 

### Need to copy this file into the same directory as the pipeline script and edit the file paths for it to run properly. 

source /apps/profiles/modules_asax.sh.dyn
module load R/4.1.0

R CMD BATCH dada2_assigntax_NBC.R

############################################
#     Bioinformatic pipeline               #
#     ITS amplicon sequences               #
#            MAPPING                       #
############################################

## Load Modules: ##
source /apps/profiles/modules_asax.sh.dyn
module load fastx/0.0.14
module load anaconda/3-2023.03
module load vsearch

## Create Directory: ##

# It creates a directory named "mapping" using the mkdir command.
mkdir ./mapping
mkdir ./otu_table
## Define Paths: ##

# It sets the paths for a samples list file (samples.txt) and a Python script (replacefastaheaders_filename.py) using variables PATH_SAMPLES_LIST and PATH_PYTHON.
# May need to change depending on where you put things. Though the python scripts directory shouldn't need to change since this is the path to the shared files. It is an absolute path.

## Loop Over Samples: #

# It then enters a loop that iterates over each sample listed in the samples.txt file.

## Inside the loop: ##

# It prints the current sample being processed.
# Uses the fastq_to_fasta command to convert a FASTQ file (merged/${sample}_merged.fastq) to a FASTA file (mapping/${sample}_merged.fasta), and it adds verbosity (-v) and
# specifies that we want to take everything including ambiguous characters (-n).
# Calls a Python script (replace_fastaheaders2.py) on the generated FASTA file from the shared directory.
# The python script replaces the fasta headers with the file name

for sample in $(cat $SAMPLES_LIST)
do

echo "On sample: $sample"
    fastq_to_fasta -i trimmed/${sample}_trimmed.fastq -o mapping/${sample}_trimmed.fasta -v -n

    # have to replace the beginning of the fasta headers with the file name for mapping. Otherwise we # get one sample with all the read counts, which is not what we want.
    python ~/noel_shared/python_scripts/replace_fastaheaders2.py mapping/${sample}_trimmed.fasta --output-dir mapping/

done

## Concatenate FASTA Files: ##

# After processing all samples, it concatenates all files with the suffix _newheaders.fasta into a single file named mapping/demultiplexed_new.fasta using the cat command.

cat mapping/*_newheaders.fasta > mapping/demultiplexed_new.fasta

## Align Demultiplexed Reads: ##

# Uses the vsearch command to perform a global alignment of the demultiplexed reads (mapping/demultiplexed_new.fasta) against a reference database (clustered/otus.fasta).
# The output is written to a file named otu_table_ITS, and it sets the identity threshold to 97% (-id 0.97).

vsearch -usearch_global mapping/demultiplexed_new.fasta -db clustered/otus.fasta -strand plus -id 0.97 -otutabout otu_table/otu_table_ITS_Fungi.txt


