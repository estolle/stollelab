#!/usr/bin/bash -l
############################################################################
# repeat and low complexity finding in genomic fasta
# Eckart Stolle, 10-2023
############################################################################
# usage: ~/scripts/repeats.sh $GENOMICFASTA $CPUs $SPECIESNAME

if [ $# -ne 3 ]; then
    echo $0: usage: ./repeats.sh $GENOMICFASTA $CPUs $SPECIESNAME 
	echo "\GENOMICFASTA: genome.fa, genome.fasta, genome.fa.gz"
	echo "\nCPUs: specify Threads to use"
	echo "\nSPECIESNAME: Species Name to use for the Output (no whitespace, e.g. GenusnameSpeciesname)"
    exit 1
fi

## set/get variables and store in a file to be use after activating conda env
INPUT=$1
CPUs=$2
SPECIESNAME=$3

## test if input exists
if [ -f "$INPUT" ]; then
    echo "INPUT exists"

else
    echo "INPUT does not exist: $INPUT" && exit 1
fi

## check if input is gzipped, unpack if needed, et File and Foldernames
if [[ $INPUT == *.gz ]]
then 
 echo "gzipped input, unpacking"
 pigz -dfk $INPUT
 INPUTFILE=$(echo $INPUT | rev | cut -d"/" -f1 | cut -d"." -f2- | rev) && echo $INPUTFILE
 INPUTFILENAME=$(echo -n $INPUTFILE | rev | cut -d"." -f2- | rev ) && echo $INPUTFILENAME
 INPUTFOLDER=$(echo $INPUT| rev | cut -d"/" -f2- | rev) && echo $INPUTFOLDER
else
 echo "not gzipped input"
  INPUTFILE=$(echo $INPUT | rev | cut -d"/" -f1 | cut -d"." -f1- | rev) && echo $INPUTFILE
  INPUTFILENAME=$(echo -n $INPUTFILE | rev | cut -d"." -f2- | rev ) && echo $INPUTFILENAME
  INPUTFOLDER=$(echo $INPUT | rev | cut -d"/" -f2- | rev) && echo $INPUTFOLDER
fi

CURRENTFOLDER=$PWD
cd $INPUTFOLDER

echo $INPUT > repeat-analyses.variables.txt
echo $CPUs >> repeat-analyses.variables.txt
echo $SPECIESNAME >> repeat-analyses.variables.txt
echo $CURRENTFOLDER >> repeat-analyses.variables.txt
echo $INPUTFILE >> repeat-analyses.variables.txt
echo $INPUTFILENAME >> repeat-analyses.variables.txt
echo $INPUTFOLDER >> repeat-analyses.variables.txt
echo $(date) >> repeat-analyses.variables.txt


# needs earlygrey installed via conda (and activated here)
# needs a bunch of tools: 

: '  ## start of block comment
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################

# install tools

mkdir -p ~/progz
mkdir -p ~/bin
# add ~/bin to your PATH in your .bashrc

## general sequence processing: seqtk
cd ~/progz
git clone https://github.com/lh3/seqtk.git; cd seqtk; make
ln -s $PWD/seqtk ~/bin

## tabtk
cd ~/progz
git clone https://github.com/lh3/tabtk; cd tabtk; make
ln -s $PWD/tabtk ~/bin

## bioawk
mkdir -p ~/progz && cd ~/progz
git clone https://github.com/lh3/bioawk; cd bioawk; make
ln -s $PWD/bioawk ~/bin

## trf-mod    # output columns: #ctg start end period copyNum fracMatch fracGap score entroy pattern
cd ~/progz
git clone https://github.com/lh3/TRF-mod; cd TRF-mod; make -f compile.mak
ln -s $PWD/trf-mod ~/bin

## etrf #exact tandem repeats (no mismatches)
cd ~/progz
git clone https://github.com/lh3/etrf; cd etrf; make
ln -s $PWD/etrf ~/bin

## bedtok
cd ~/progz
git clone https://github.com/lh3/bedtk; cd bedtk; make
ln -s $PWD/bedtk ~/bin

# satellite repeatfinder srf
cd ~/progz
git clone https://github.com/lh3/srf; cd srf; make
ln -s $PWD/srf ~/bin

## kmc
cd ~/progz
mkdir kmc; cd kmc
wget https://github.com/refresh-bio/KMC/releases/download/v3.2.2/KMC3.2.2.linux.x64.tar.gz
tar -xvf KMC3.2.2.linux.x64.tar.gz
ln -s $PWD/bin/kmc ~/bin
ln -s $PWD/bin/kmc_tools ~/bin
ln -s $PWD/bin/kmc_dump ~/bin

## minimap2
cd ~/progz
git clone https://github.com/lh3/minimap2; cd minimap2; make
ln -s $PWD/minimap2 ~/bin

## sdust (mdust alternative)
cd ~/progz
git clone https://github.com/lh3/sdust; cd sdust; make
ln -s $PWD/sdust ~/bin


## general sequence processing: seqkit
cd ~/progz
wget https://github.com/shenwei356/seqkit/releases/download/v2.5.1/seqkit_linux_amd64.tar.gz
tar -xf seqkit_linux_amd64.tar.gz
ln -s $PWD/seqkit ~/bin

## SSR detection with ULTRA
cd ~/progz
git clone https://github.com/TravisWheelerLab/ULTRA; cd ULTRA; cmake .; make
ln -s $PWD/ultra ~/bin

## SSR detection with TANTAN
cd ~/progz
git clone https://gitlab.com/mcfrith/tantan; cd tantan; make
ln -s $PWD/bin/tantan ~/bin

## mappability genmap
cd ~/progz
wget https://github.com/cpockrandt/genmap/releases/download/genmap-v1.3.0/genmap-1.3.0-Linux-x86_64-sse4.zip
unzip genmap-1.3.0-Linux-x86_64-sse4.zip
ln -s $PWD/genmap ~/bin

## bedtools
cd ~/progz
mkdir -p bedtools/v2.31.0; cd bedtools/v2.31.0
wget https://github.com/arq5x/bedtools2/releases/download/v2.31.0/bedtools.static
chmod 755 bedtools.static
ln -s $PWD/bedtools.static ~/bin/bedtools

## mosdepth
cd ~/progz
wget https://github.com/brentp/mosdepth/releases/download/v0.3.5/mosdepth
chmod 755 mosdepth
ln -s $PWD/mosdepth ~/bin

## bwa-mem2
cd ~/progz
curl -L https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.2.1/bwa-mem2-2.2.1_x64-linux.tar.bz2 \
  | tar jxf -
cd bwa-mem2-2.2.1_x64-linux
ln -s $PWD/bwa-mem2.avx2 ~/bin/bwa
ln -s $PWD/bwa-mem2.avx2 ~/bin/

#genometools
cd ~/progz
wget http://genometools.org/pub/genometools-1.6.2.tar.gz
tar -xf genometools-1.6.2.tar.gz
cd genometools-1.6.2
make threads=yes 64bit=yes
# if error cairo.h missing in that case try:
# make threads=yes 64bit=yes cairo=no
# or install systemwide cairo (root required) sudo apt-get install libcairo2-dev libpango1.0-dev
ln -s $PWD/bin/gt ~/bin
ln -s $PWD/bin/lua ~/bin

###### EarlGrey installation via conda
## in the earlygrey exec change/check the cd-hit-est command (see github issue)
## this sinstallation was done in a systemwide anaconda via sudo priviledge, but adaptation 
## for miniconda should be possible
## in this installation I could not run the ./configure in RepeatMasker inside this conda to be able to integrate
## a RepBase repeat Database. I installed RM separately, ran the configuration there so that the library was generated
## and then I copied it here into the Libraries folders
## mamba was installed already

## activate conda
cd /opt/anaconda3-2023.07-2-Linux-x86_64/
su anaconda
export PATH="/opt/anaconda3-2023.07-2-Linux-x86_64/bin:$PATH"
source /opt/anaconda3-2023.07-2-Linux-x86_64/bin/activate

## create conda env in custom location (i.e. not in home where space is limited)
conda create --prefix=/scratch/progz/conda_envs/earlgrey
conda info --envs  ##check all your envs
conda clean --all --force-pkgs-dirs  ##if you need to re-try the install due to error

## activate env and install packages
conda activate /scratch/progz/conda_envs/earlgrey
mamba install -c conda-forge -c bioconda earlgrey

# find RepeatMasker location within the env
which RepeatMasker
ll /scratch/progz/conda_envs/earlgrey/bin/RepeatMasker
cd /scratch/progz/conda_envs/earlgrey/share/RepeatMasker/RepeatMasker
mv Libraries Libraries.bak
#copy Libraries folder from a fully configured RepeatMasker install with full library
# modify this source PATH to your own RM installation+configured libs
cp -r /scratch/progz/RepeatMasker/Libraries .


## install repeatmasker separately, configure and include the RepBase library
# this can be then used to copy over (above)
# before running ./configure, softlink any DFAM.hmm or RepBaseRMedition.embl


######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
'  ## end of block comment

## activate conda env 8needs additional steps to enable conda in bash scripts)
export PATH="/opt/anaconda3-2023.07-2-Linux-x86_64/bin:$PATH" 
eval "$(conda shell.bash hook)"
#source /opt/anaconda3-2023.07-2-Linux-x86_64/etc/profile.d/conda.sh
#export PATH="/opt/anaconda3-2023.07-2-Linux-x86_64/bin:$PATH"
#source /opt/anaconda3-2023.07-2-Linux-x86_64/bin/activate
conda activate earlgrey


## reasssign variables (after conda activation)

INPUT=$(cat repeat-analyses.variables.txt | sed -n '1p')
CPUs=$(cat repeat-analyses.variables.txt | sed -n '2p')
SPECIESNAME=$(cat repeat-analyses.variables.txt | sed '3q;d')
CURRENTFOLDER=$(cat repeat-analyses.variables.txt | sed '4q;d')
INPUTFILE=$(cat repeat-analyses.variables.txt | sed '5q;d')
INPUTFILENAME=$(cat repeat-analyses.variables.txt | sed '6q;d')
INPUTFOLDER=$(cat repeat-analyses.variables.txt | sed '7q;d')
STARTDATE=$(cat repeat-analyses.variables.txt | sed '8q;d')


## make outputfolder
SPECIES="$SPECIESNAME"
OUTPUT="$INPUTFILE.earlgrey"
REF="$INPUTFILE"
echo "$OUTPUT" >> repeat-analyses.variables.txt


## run EarlGrey pipeline

## earlygrey
time earlGrey -g $REF -s $SPECIES -o $OUTPUT -t $CPUs -c yes

## test if output exists
OUTPUTFILE=$(ls $OUTPUT/$SPECIES'_EarlGrey'/$SPECIES'_summaryFiles'/$SPECIES.filteredRepeats.bed)

if [ -s "$OUTPUTFILE" ]; then
    echo "OUTPUTFILE exists, continuing"
else
    echo "earlgrey output does not exist or is empty, restarting"
    # exit 1
    #mv TS* folder inside the strainer folder and restart earlygrey to rerun cdhit
    mv $OUTPUT/$SPECIES'_EarlGrey'/$SPECIES'_strainer'/'TS_'$SPECIES'-families'.* $OUTPUT/$SPECIES'_EarlGrey'/
    time earlGrey -g $REF -s $SPECIES -o $OUTPUT -t $CPUs -c yes
fi


## process some files
cat $OUTPUT/$SPECIES'_EarlGrey'/$SPECIES'EarlGrey'.log | grep -v "% completed" | grep -v "\[0m" | grep -v "########" | grep -v 'running on ctg' > $OUTPUT/$SPECIES'_EarlGrey'/$SPECIES'EarlGrey'.short.log
cat $OUTPUT/$SPECIES'_EarlGrey'/$SPECIES'_summaryFiles'/$SPECIES.filteredRepeats.bed | bioawk -t '{print $0,$3-$2}' > $OUTPUT/$SPECIES'_EarlGrey'/$SPECIES'_summaryFiles'/$SPECIES.filteredRepeats.withSize.bed
cat $OUTPUT/$SPECIES'_EarlGrey'/$SPECIES'_summaryFiles'/$SPECIES.filteredRepeats.withSize.bed | bioawk '{sum+=$7;} END{print sum;}' > $OUTPUT/$SPECIES'_EarlGrey'/$SPECIES'_summaryFiles'/$SPECIES.filteredRepeats.withSize.bed.sum


## cleanup
if [ -s "$OUTPUTFILE" ]; then
echo "EarlyGrey Output exists, cleaning up"
cd $OUTPUT/$SPECIES'_EarlGrey'/$SPECIES'_Database'
rm -f $SPECIES.n*
pigz *.stk *.fa
cd ../$SPECIES'_RepeatMasker_Against_Custom_Library' 
rm -f $REF.prep.align $REF.prep.cat.gz
pigz *.out *.masked
cd ../$SPECIES'_mergedRepeats'/looseMerge
pigz *.gff* *.gff
cd ../../$SPECIES'_RepeatModeler'/RM_*/
pigz *.stk
rm -rf ./round-* 
cd $INPUTFOLDER
else
echo "EarlyGrey OUTPUTFILE doesnt exist, perhaps the tool failed, exiting"
# && exit 1
 fi


# copy summary file(s))
cp $OUTPUT/$SPECIES'_EarlGrey'/$SPECIES'_summaryFiles'/$SPECIES.filteredRepeats.withSize.bed ./$REF.earlgrey.filteredRepeats.withSize.bed
cp $OUTPUT/$SPECIES'_EarlGrey'/$SPECIES'_summaryFiles'/$SPECIES.filteredRepeats.withSize.bed.sum ./$REF.earlgrey.filteredRepeats.withSize.bed.sum
cat $REF.earlgrey.filteredRepeats.withSize.bed | grep -vP 'Simple_repeat|Low_complexity|Satellite|Unknown' > $REF.earlgrey.filteredRepeats.withSize.TEonly.bed 
cat $REF.earlgrey.filteredRepeats.withSize.TEonly.bed | bioawk '{sum+=$7;} END{print sum;}' > $REF.earlgrey.filteredRepeats.withSize.TEonly.bed.sum

printf "EarlGrey\t" > repeats.stats.txt
cat $REF.earlgrey.filteredRepeats.withSize.bed.sum >> repeats.stats.txt
printf "EarlGreyTEonly\t" >> repeats.stats.txt
cat $REF.earlgrey.filteredRepeats.withSize.TEonly.bed.sum >> repeats.stats.txt

# Fasta Index
samtools faidx $REF
cat $REF.fai | cut -f1,2 > $REF.length

# occurences of N
seqtk cutN -g -n 1 $REF | awk '{print $0,$3-$2}' > $REF.NNN.bed
cat $REF.NNN.bed | bioawk '{sum+=$4;} END{print sum;}' > $REF.NNN.bed.sum
printf "N\t" >> repeats.stats.txt
cat $REF.NNN.bed.sum >> repeats.stats.txt

# long homopolymer runs
seqtk hrun $REF | awk '{print $0,$3-$2}' > $REF.hrun.bed
cat $REF.hrun.bed | bioawk '{sum+=$5;} END{print sum;}' > $REF.hrun.bed.sum
printf "hrun\t" >> repeats.stats.txt
cat $REF.hrun.bed.sum >> repeats.stats.txt

# sdust
sdust -w 64 -t 20 $REF | awk '{print $0,$3-$2}' > $REF.sdust.bed
cat $REF.sdust.bed | bioawk '{sum+=$4;} END{print sum;}' > $REF.sdust.bed.sum
printf "sdust\t" >> repeats.stats.txt
cat $REF.sdust.bed.sum >> repeats.stats.txt

# mappability
# FYI mappability=1: k-mer occurs once with up to e errors, low mappability: k-mer belongs to a repetitive region
genmap index --fasta-file $REF --index $REF.genmap.index
genmap map -K 35 -E 2 --threads $CPUs --bedgraph --wig --txt --output $REF.genmap -I $REF.genmap.index
#ls $REF.genmap*
cat $REF.genmap.bedgraph | awk '$4<=0.25' | cut -f1-3 | bedtools merge -d 10 | bioawk -t '{print $0,$3-$2}' > $REF.genmap.0.25.bed
cat $REF.genmap.0.25.bed | bioawk '{sum+=$4;} END{print sum;}' > $REF.genmap.0.25.bed.sum
cat $REF.genmap.bedgraph | awk '$4<=0.2' | cut -f1-3 | bedtools merge -d 10 | bioawk -t '{print $0,$3-$2}' > $REF.genmap.0.2.bed
cat $REF.genmap.0.2.bed | bioawk '{sum+=$4;} END{print sum;}' > $REF.genmap.0.2.bed.sum
pigz $REF.genmap.bedgraph $REF.genmap.wig $REF.genmap.txt
printf "genmap.0.2\t" >> repeats.stats.txt
cat $REF.genmap.0.2.bed.sum >> repeats.stats.txt
rm -rf $REF.genmap.index


## SSRs
tantan -f4 $REF | bioawk -t '{print $0,$3-$2}' > $REF.SSRs.tantan.bed
cat $REF.SSRs.tantan.bed | bioawk '{sum+=$8;} END{print sum;}' > $REF.SSRs.tantan.bed.sum
printf "SSRs.tantan\t" >> repeats.stats.txt
cat $REF.SSRs.tantan.bed.sum >> repeats.stats.txt

trf-mod -p 300 $REF | bioawk -t '{print $0,$3-$2}' > $REF.SSRs.trfmod.bed
cat $REF.SSRs.trfmod.bed | bioawk '{sum+=$4;} END{print sum;}' > $REF.SSRs.trfmod.bed.sum
printf "SSRs.trfmod\t" >> repeats.stats.txt
cat $REF.SSRs.trfmod.bed.sum >> repeats.stats.txt

etrf -m 200 -l 13 $REF | bioawk -t '{print $0,$3-$2}' > $REF.SSRs.etrf.bed
cat $REF.SSRs.etrf.bed | bioawk '{sum+=$7;} END{print sum;}' > $REF.SSRs.etrf.bed.sum
printf "SSRs.etrf\t" >> repeats.stats.txt
cat $REF.SSRs.etrf.bed.sum >> repeats.stats.txt

ultra --minunit 4 --minlen 10 --period 300 --threads $CPUs $REF | grep -v "^\"" | grep -v "^\{" | grep -v "^\}" > $REF.SSRs.ultra.bed
cat $REF.SSRs.ultra.bed | bioawk '{sum+=$9;} END{print sum;}' > $REF.SSRs.ultra.bed.sum
cat $REF.SSRs.ultra.bed.sum
printf "SSRs.ultra\t" >> repeats.stats.txt
cat $REF.SSRs.ultra.bed.sum >> repeats.stats.txt

#srf
mkdir -p srf_tmp_dir
kmc -fm -k171 -t$CPUs -ci20 -cs100000 $REF $REF.count.kmc ./srf_tmp_dir
kmc_dump $REF.count.kmc $REF.count.txt
srf -p $REF $REF.count.txt > $REF.srf.fa
# enlong short contigs for mapping
srfutils.js enlong $REF.srf.fa > $REF.srf.enlong.fa
minimap2 -c -N1000000 -f1000 -r100,100 -t$CPUs $REF.srf.enlong.fa $REF > $REF.srf.paf
# generate non-redundant mapping regions
srfutils.js paf2bed $REF.srf.paf > $REF.srf.bed
srfutils.js bed2abun $REF.srf.bed > $REF.srf.bed.abundance
trf-mod $REF.srf.fa > $REF.srf.decomposed.fa
# TRF-mod: decompose HORs to monomers; #many similar seq are mapped in same repeat? rerun and filter customly: minimap2 -c -N1000 <(./srfutils.js enlong -d srf.fa) srf.fa
rm -rf srf.tmp_dir
cat $REF.srf.bed | bioawk -t '{print $0,$3-$2}' > $REF.srf.with.size.bed
cat $REF.srf.with.size.bed | bioawk '{sum+=$9;} END{print sum;}' > $REF.srf.with.size.bed.sum
printf "Sats.srf\t" >> repeats.stats.txt
cat $REF.srf.with.size.bed.sum >> repeats.stats.txt


cat $REF.earlgrey.filteredRepeats.withSize.TEonly.bed \
$REF.NNN.bed \
$REF.sdust.bed \
$REF.SSRs.ultra.bed \
$REF.SSRs.etrf.bed \
$REF.SSRs.trfmod.bed \
$REF.SSRs.tantan.bed \
$REF.srf.with.size.bed \
$REF.genmap.0.2.bed \
$REF.hrun.bed \
  | cut -f1-3 | sort -k1,1 -k2,2n \
  | bedtools slop -i /dev/stdin -g $REF.length -b 3 \
  | bedtools merge -d 6 | bedtools sort \
  | bioawk -t '{print $0,$3-$2}' \
  | bgzip -f -@ $CPUs -c /dev/stdin > $REF.repetitive.bed.gz
zcat $REF.repetitive.bed.gz | bioawk '{sum+=$4;} END{print sum;}' > $REF.repetitive.bed.sum
printf "Merged.repetitive\t" >> repeats.stats.txt
cat $REF.repetitive.bed.sum >> repeats.stats.txt

cat repeats.stats.txt

ENDTIME=$(date)
echo "$ENDTIME" >> repeat-analyses.variables.txt

echo "finished repeat analysis, merged bed file is here: \n $REF.repetitive.bed"
echo "started: "$STARTTIME
echo "finished at: "$ENDTIME
cd $CURRENTFOLDER

#conda deactivate
echo "done"

