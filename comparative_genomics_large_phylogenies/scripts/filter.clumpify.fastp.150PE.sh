#!/bin/bash
############################################################################
# readfiltering with clumpify and fastp
# Eckart Stolle, 2023
############################################################################
# usage: ~/scripts/filter.clumpify.fastp.150PE.sh $INPUTFWD $CPUs

##RNAseq data
#cd /scratch/ek/reads/2021.Camptopoeum.friesei/ERGA.RNAseq.2022
CPUs=16
#INPUTFWD="/scratch/ek/reads/2021.Camptopoeum.friesei/ERGA.RNAseq.2022/DE24NGSUKBR128751_Cf2_S3_R1_001.fastq.gz"
#/scratch/ek/reads/2021.Camptopoeum.friesei/ERGA.RNAseq.2022/DE24NGSUKBR128751_Cf2_S3_R2_001.fastq.gz
#INPUTFWD="/scratch/ek/reads/2021.Camptopoeum.friesei/ERGA.RNAseq.2022/DE94NGSUKBR128752_Cf3_S4_R1_001.fastq.gz"
#/scratch/ek/reads/2021.Camptopoeum.friesei/ERGA.RNAseq.2022/DE94NGSUKBR128752_Cf3_S4_R2_001.fastq.gz


if [ $# -ne 2 ]; then
    echo $0: usage: ./filter.clumpify.fastp.150PE.sh INPUTFWD CPUs 
	echo "\nINPUT: fwd reads, format PATH/nnnn_R1_001.fastq.gz, rev reads would be PATH/nnnn_R2_001.fastq.gz"
	echo "\nCPUs: specify Threads to use, note: fastp uses max 16 Threads"
    exit 1
fi

CPUs=16
AVG_QUAL_THRESHOLD=30
END_QUAL_THRESHOLD=20
L=80
OPTICALDIST=12000
##dupedist instead dist HiSeq2000:40 HiSeq4000:2500 Novaseq:12000

## specific adapter sequences for removal
#FWD_ADAPTERS="$HOME/scripts/illumina.fwd.adapters.fa"
#REV_ADAPTERS="$HOME/scripts/illumina.rev.adapters.fa"



#set/get variables
INPUTFWD=$1
CPUs=$2

INPUTFILE=$(echo $INPUTFWD | rev | cut -d"." -f 2- | cut -d"_" -f 3- | cut -d"/" -f 1 | rev)
echo $INPUTFILE
INPUTFOLDER=$(echo $INPUTFWD | rev | cut -d"/" -f 2- | rev)
echo $INPUTFOLDER


FWD="$INPUTFOLDER/$INPUTFILE""_R1_001.fastq.gz"
REV="$INPUTFOLDER/$INPUTFILE""_R2_001.fastq.gz"
ls $FWD $REV

## test if these exist
if [ -f "$FWD" ]; then
    echo "INPUT FWD fastq exists"
else
    echo "INPUT FWD fastq exists: $INPUTFWD" && exit 1
fi 

echo "filtering reads with clumpify"
clumpify.sh in1=$FWD in2=$REV out1=$INPUTFOLDER/$INPUTFILE.optdepup.R1.fastq.gz out2=$INPUTFOLDER/$INPUTFILE.optdepup.R2.fastq.gz dedupe optical dupedist=$OPTICALDIST -Xmx30g passes=1 subs=1 k=31 spantiles=f 2> >(tee $INPUTFOLDER/$INPUTFILE.optical_duplicates.log >&1)
echo "clumpify_done"

echo "running readfiltering with fastp"
fastp --thread $CPUs --length_required $L --average_qual $AVG_QUAL_THRESHOLD --cut_right --cut_window_size 10 --cut_mean_quality $END_QUAL_THRESHOLD --n_base_limit 1 -z 7 --in1 $INPUTFOLDER/$INPUTFILE.optdepup.R1.fastq.gz --in2 $INPUTFOLDER/$INPUTFILE.optdepup.R2.fastq.gz --out1 $INPUTFOLDER/$INPUTFILE.filtered.R1.fq.gz --out2 $INPUTFOLDER/$INPUTFILE.filtered.R2.fq.gz --json $INPUTFOLDER/$INPUTFILE.filtered.fastp.report.json --html $INPUTFOLDER/$INPUTFILE.filtered.fastp.report.html --report_title $INPUTFILE
echo "fastp_done"







## old code when running folders

#ls -1 $INPUTFOLDER/*_[1,2].fastq.gz | rev | cut -d"." -f 2- | cut -d"_" -f 2- | cut -d"/" -f 1 | rev | sort | uniq > $OUTPUTFOLDER/samples.lst
#ls -1 $INPUTFOLDER/*_[1,2].fastq.gz | rev | cut -d"." -f 2- | cut -d"_" -f 2- | cut -d"/" -f 1 | rev | sort | uniq > $OUTPUTFOLDER/samples.lst
#LISTtoPROCESS="$OUTPUTFOLDER/samples.lst"
#N=$(wc -l $LISTtoPROCESS | cut -d" " -f1)
#echo "processing $N samples (read-pair files)"

#ls -1 *.fastq.gz | cut -d "_" -f 1,2 | uniq > samples.lst
#cat samples.lst | parallel "echo {}; ls {}*R1*.fastq.gz"
#cat samples.lst | parallel "echo {}; mv {}*R1*.fastq.gz {}.R1.fastq.gz"
#cat samples.lst | parallel "echo {}; mv {}*R2*.fastq.gz {}.R2.fastq.gz"
##readfilter
#cat samples.lst | parallel "echo {}; ls {}.R1.fastq.gz {}.R2.fastq.gz"
#cat samples.lst | parallel -k -j5 "echo {}; fastp --thread $CPUs --length_required 50 --average_qual 30 --cut_right --cut_window_size 10 --cut_mean_quality 25 --n_base_limit 1 -z 7 --in1 {}.R1.fastq.gz --in2 {}.R2.fastq.gz --out1 {}.filtered.R1.fq.gz --out2 {}.filtered.R2.fq.gz --json {}.filtered.fastp.report.json --html {}.filtered.fastp.report.html --report_title {}"





