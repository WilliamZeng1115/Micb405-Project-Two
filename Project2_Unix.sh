
################################################
#!/bin/bash

clear
echo `date`
echo "MICB 405 2018W1 Project 2 - Group 4"
echo "Run Prokka:"

cd ~
mkdir project2
cd ~/project2

mkdir PROKKA

#Global Variables
DEPTH=10 


ANNOTE_PATH=/projects/micb405/resources/project_2/2018/SaanichInlet_10m/MetaBAT2_SaanichInlet_10m/gtdbtk_output
MAGS_PATH=/projects/micb405/resources/project_2/2018/SaanichInlet\_10m/MetaBAT2\_SaanichInlet_10m/MedQPlus_MAGs

#awk {print $2} awk print the second column of the input
#awk -F";" '{ print $1 }' awk -F specify the file separator, thus this takes everything before ';'
#sed 's/d__//g' replace "d__" with empty character

for f in $MAGS_PATH/*.fa
do 
bin=${f#*.};
bin=${bin%.*};
tax=$(grep -w $bin $ANNOTE_PATH/gtdbtk.*.classification_pplacer.tsv | awk '{ print $2 }' | awk -F";" '{ print $1 }'| sed 's/d__//g');
prokka --kingdom $tax --outdir ~/project2/PROKKA/${bin}/ --force --prefix ${bin}_SaanichInlet_MAG_ORFs $f
done

>SannichInlet_MAG_ORFs.faa
for f in ~/project2/PROKKA/*/*.faa
do
echo $f
cat $f >>SannichInlet_MAG_ORFs.faa
done


#https://www.genome.jp/kaas-bin/kaas_main?mode=queryinfo&id=1542469991&key=H2TrLKPk
#ID :	1542469991
#query name :	query
#program :	GHOSTX
#method :	SBH
#GENES data set :	hsa, dme, cel, ath, sce, cho, eco, nme, hpy, rpr
#bsu, lla, cac, mge, mtu, ctr, bbu, syn, bth, dra
#aae, mja, ape

######################################
echo "Transcriptional Activity, RPKM"

TXN_PATH=/projects/micb405/resources/project_2/2018/Metatranscriptomes
DEPTH=10m
mags=(155 192 33 235 120 65 84 262 247 293 378)


cd ~/project2

mkdir SAM_PROKKA
mkdir MAG_RPKM

for i in ${mags[@]}
do
	f=~/project2/PROKKA/$i/$i_*.ffn
	mkdir ~/project2/MAG_RPKM/$i
	bwa index $f
	for q in $TXN_PATH/*$DEPTH*.fastq*
	do
		sid=${q##*/}
		sid=${sid%%_*} 
		bwa mem -t 2 $f $q > ~/project2/SAM_PROKKA/${sid}_${i}_MAG_ORFs.sam
		/projects/micb405/resources/project_2/2018/rpkm -c $f -a ~/project2/SAM_PROKKA/${sid}_${i}_MAG_ORFs.sam -o ~/project2/MAG_RPKM/${i}/${i}_${sid}_MAG_ORFs_RPKM.csv

	done
done


