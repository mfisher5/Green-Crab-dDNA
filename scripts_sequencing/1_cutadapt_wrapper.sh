#!/bin/bash

#######
# Purpose: Run the Cutadapt wrapper script, which calls R/run.cutadapt.sh
# Command Line: bash scripts/1_cutadapt_wrapper.sh {folder with fastqs} {metadata_file} {output_folder}
# Environment: git bash on Windows 10, GNU bash version  4.4.23(2)-release (x86_64-pc-msys)
# Author: M Fisher via Eily Allen
# Date Written: 2022-02-02
# Last Updated: 2022-06-16
#######

MAIN_DIR="$(dirname "$0")"
pwd
#mkdir "${OUTPUT_DIR}"
cat ${1}

# set params
METADATA=($(awk -F',' -v COLNUM=2 \
  'NR>1 {  print $COLNUM }' ${1} \
   ))
MINLENGTH=($(awk -F',' -v COLNUM=4 \
  'NR>1 {  print $COLNUM }' ${1} \ 
  ))

#Capture one value of the params file
FASTQFOLDER=($(awk -F',' -v COLNUM=1 \
  'NR>1 {  print $COLNUM }' ${1} \
   ))
echo "fastq files will be read in from this folder:"
echo  "${FASTQFOLDER}"
echo "---"

OUTPUTFOLDER=($(awk -F',' -v COLNUM=3 \
  'NR>1 {  print $COLNUM }' ${1} \
   ))
echo "and trimmed fastq files will be saved into this folder:"
echo  "${OUTPUTFOLDER}"
echo "---"
echo "---"
echo "now running cutadapt script... output will be saved to a log file."

LOGFILE=($(date +%y-%b-%d_cutadapt-log.txt))
LOGDIR="data/cutadapt/"
LOGDIR+=${LOGFILE}



bash R/run.cutadapt.sh "${FASTQFOLDER}" "${METADATA}" "${OUTPUTFOLDER}" "${MINLENGTH}" >> "${LOGDIR}"