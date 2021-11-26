#!/bin/bash

usage()
{
	cat << EOF
	usage: $0 optionsi
	Ensure you are in chos 7 before running the command
	
	This command takes in a pacbio bam file and a pacbio assmebly to polish it with arrow
	It will run arrow multiple times to iteratively polish it
	This will stop after a user set amount of times or when no changes were carried out by the last polishing
	OPTIONS:
		-h		Show this
		-s <prefix>	sample_name, this will be used for directory and file names
		-f <FILE>	Assembly fasta file
		-b <FILE>	Pacbio subread xml file
		-i <INT>	Number of max polishing iterations
		-t <INT>	Number of threads/workers for pbalign and arrow commands 
		-o <DIR>	output directory, ensure this is created before running
EOF
}

#Variables
SAMPLE=
FASTA=
BAM=
MAX_ITERATIONS=
NUM_THREADS=
OUTPUT_DIR=

#Get options
while getopts "hxs:f:b:i:t:o:" OPTION
do
	case $OPTION in 
		#help
		h)
			usage
			exit 1
			;;
		#Prefix of files
		s)
			SAMPLE=$OPTARG
			;;
		#Fasta assembly file
		f)
			FASTA=$OPTARG
			;;
		#BAM subread file
		b)
			BAM=$OPTARG
			;;
		#Max number of polishign iterations
		i)
			MAX_ITERATIONS=$OPTARG
			;;
		#Num of threads for commands to use
		t)
			NUM_THREADS=$OPTARG
			;;
		#OUTPUT directory
		o)
			OUTPUT_DIR=$OPTARG
			;;
		
	esac
done

FASTA=$(readlink -f $FASTA)
BAM=$(readlink -f $BAM)
OUTPUT_DIR=$(readlink -f $OUTPUT_DIR)

#Check if files and direcotires exist
echo "###This Genome polishing is for sample:" ${SAMPLE}
echo "###Check Fasta assmebly"
ls $FASTA
echo "###Check Bam subread"
ls $BAM
MAIN_DIR=${OUTPUT_DIR}/${SAMPLE}/
#Make subdirectory for specific sample
mkdir ${MAIN_DIR}
cd ${MAIN_DIR}
echo "##Check Output Directory##"
ls $MAIN_DIR


#Number the iteration is at
iteration=1
#Set temp assmebly variable that will be changed
ASSEMBLY=${FASTA}
#Below is a temp variable which will allow loop to start but may change in the loop to stop loop
#This will occur if the last polishing iteration found no errors
ERRORS_FOUND="1"

cd ${MAIN_DIR}
#Loop for arrow polishing
while [ "$iteration" -le "${MAX_ITERATIONS}" ]
do
        #Set path of output pbalign.bam and arrow outputs
        PBMM_BAM=${MAIN_DIR}/${SAMPLE}_iteration_${iteration}_pbmm.bam
        ARROW_FASTA=${MAIN_DIR}/${SAMPLE}_iteration_${iteration}_arrow_consensus.fasta
        ARROW_FASTQ=${MAIN_DIR}/${SAMPLE}_iteration_${iteration}_arrow_consensus.fastq
	#Index Assembly for pbmm2 alignment
	~pbio7/smrtlink/smrtcmds/bin/pbmm2 index ${ASSEMBLY} ${ASSEMBLY}.mmi
	#pbmm2 alignment
	~pbio7/smrtlink/smrtcmds/bin/pbmm2 align --sort --preset "SUBREAD" ${ASSEMBLY}.mmi ${BAM} ${PBMM_BAM}
        #Create fasta companion file
        samtools faidx ${ASSEMBLY}
        #Index pbalign bam file
        ~pbio7/smrtlink/smrtcmds/bin/pbindex ${PBMM_BAM}
        #Arrow polishing
        ~pbio7/smrtlink/smrtcmds/bin/quiver --referenceFilename ${ASSEMBLY} -o ${ARROW_FASTA} -o ${ARROW_FASTQ} -j ${NUM_THREADS} ${PBMM_BAM}
        #Assign assembly variable to new assmebly for next iteration
        ASSEMBLY=$ARROW_FASTA
        #Find out how many variants were found
        
        #This loop does not appear to work so i will annotate it out for the moment
	#if ["${ERRORS_FOUND}" -eq 0]; then
        #        break
        #fi
        #Add one to counter
        iteration=$(( $iteration + 1 ))
done

#remove temp directories
rm -r ${MAIN_DIR}/aligned_split_bams_iter_*

