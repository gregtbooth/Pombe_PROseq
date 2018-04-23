#source  ~/.bash_profile
#====================
#PRO-seq shell script
#Modified for simultaneous mapping of spiked in organism data
#Ver. 2.0	     
#====================
#==========================================================================================================================================
#Ingredients					      
#1. Your FASTA, rDNA (ebwt), genome (ebwt) *ebwt can be generated by doing "bowtie-build infile outfile"*  				      
#2. One python script: passfilter.py
#3. fastx_tools					      
#4. bowtie					      
#5. Samtools
#6. BEDTools				              
#===========================================================================================================================================
#Usage
#1. Input file path of necessary files and scripts
#2. Run the script by "sh SCRIPTPATH"
#NOTE: If error arises in quality filter, change -o to -q 
#(optional) if no rDNA, 
#hashtag or remove all the lines in between "++++++++++++++" 
#AND hashtage or remove one "done" in the end of the script
#===========================================================================================================================================
#Input
## single file that must be parsed using fastx-toolkit to give individual ePROseq fastq files (sample names provided in barcode file)
ePROseqfastq="/home/macproadmin/users/Greg/FisherCollaboration/11-28-17_Sequencing/fastq/9031_7157_66639_H3KH2BGX5_FisherPoolA_1-5_ATCACG_R1.fastq.gz"

SpikeInrDNA="/home/macproadmin/users/GROseq_Utils/refGenomes/Cerevisiae_S288C/rDNA/rDNA"

rDNA="/home/macproadmin/users/GROseq_Utils/refGenomes/pombe/rRNA/Pombe_rRNA.fa"

genome="/home/macproadmin/users/GROseq_Utils/refGenomes/CombinedGenomes/Cerevisiae_Pombe/Combined_SC_SP_genomes"

Combinedchrominfo="/home/macproadmin/users/GROseq_Utils/annotations/CombinedGenomes/Cerevisiae_Pombe/SC_SPcombinedChrSizes.txt"

SCchrominfo="/home/macproadmin/users/GROseq_Utils/annotations/sacCer3/SacCerChr.sizes_forPipeline.txt"

SPchrominfo="/home/macproadmin/users/GROseq_Utils/annotations/pombe/PombeChrSizes_forPipeline.txt"
#chrominfo format example = chr1 52949123 (with no title)

bgToBigWig="/home/macproadmin/users/GROseq_Utils/mkbigWig/bedGraphToBigWig"

nmer="43" # 36 + 7nt barcode at beginning

barcodes="/home/macproadmin/users/Greg/FisherCollaboration/11-28-17_Sequencing/pipeline/Barcodes_samples1-5.txt"

#===========================================================================================================================================
gunzip $ePROseqfastq

cat ${ePROseqfastq%%.*}.fastq | fastx_barcode_splitter.pl --bcfile ${barcodes} --bol --mismatches 2 --partial 2 --prefix ePROseq_ --suffix ".fastq"
	## barcode file is a text file that gives the barcode sequence and the name of the files, --bol matches barcodes at the beginning of sequences, 2 mismatches allowed and also allows 2 partial overlap of barcodes and prefix and suffix are added to the name of the individual barcodes in barcode files.
gzip ${ePROseqfastq%%.*}.fastq

# after splitting files above, will have new new fastq files with names provided in the barcode_12 file
sampleList="
ePROseq_Pombe_Fisher_WT_rep1
ePROseq_Pombe_Fisher_dis2Del_rep1
ePROseq_Pombe_Fisher_dis2_11_rep1
ePROseq_Pombe_Fisher_dis2_T316A_rep1
ePROseq_Pombe_Fisher_dis2_T316D_rep1
"

for file in ${sampleList}
do

echo ${file} clipping...

#Takes _q.temp file and clips any residual 3'adapter sequence which could be all or part of TGGAATTCTCGGGTGCCAAGG and generated anothe temporary file called _c.tmp
fastx_clipper -i ${file%%.*}.fastq -o _c.tmp -a TGGAATTCTCGGGTGCCAAGG -l 15   

cat _c.tmp | awk '{if(NR%4==2) print length($1)}' > ${file%%.*}_clipLength.txt
cat ${file%%.*}_clipLength.txt | sort -n | uniq -c | awk '{print $2 " " $1}' > ${file%%.*}_clipLengthCounts.txt
rm ${file%%.*}_clipLength.txt

#Takes _c.tmp file and trims the 7 nt barcode off of sequences and then sets max length to nmer number of bases.
echo trimming first 7 nt from reads... inline barcode
fastx_trimmer -f 8 -l ${nmer} -i _c.tmp -o _c1.tmp

rm _c.tmp

echo ${file} reverse complementing...
#Takes the _c1.tmp file and flips the read because the illumnia sequences the DNA from 5'end so the adapters in PRO-seq is reversed.
fastx_reverse_complement -i _c1.tmp -o _3.tmp    

rm _c1.tmp

#+++++++++++++++++++++++++++++++++++++++++++++++++++++HASHTAG IF NO rDNA+++++++++++++++++++++++++++++++++++++++++++++++++++++++++

for rib in ${SpikeInrDNA}
do

echo ${file%%.*} aligning to SpikeInribosomal DNA...

#Generate 3 files: ribUnalign (not aligned to rDNA sequences), ribAlign (aligned to rDNA), _rib.log (bowtie stats)
#bowtie -p3 -v2 -M1 -q --sam --un ${file}_ribUnalign ${rib} _3.tmp ${file}_ribAlign 2>&1| tee ${file}_rib.log   
bowtie -p3 -k2 -m1 -q --sam --un ${file%%.*}_SpikeInribUnalign ${rib} _3.tmp ${file%%.*}_SpikeInribAlign 2>&1| tee ${file%%.*}SpikeIn_rib.log      
#cp -avr ${file}_ribAlign /home/macproadmin/users/Greg/Ribosom_alignFiles/
rm ${file%%.*}_SpikeInribAlign

#+++++++++++++++++++++++++++++++++++++++++++++++++++++HASHTAG IF NO rDNA+++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 


for rib in ${rDNA}
do

echo ${file} aligning to ribosomal DNA...

#Generate 3 files: ribUnalign (not aligned to rDNA sequences), ribAlign (aligned to rDNA), _rib.log (bowtie stats)
#bowtie -p3 -v2 -M1 -q --sam --un ${file}_ribUnalign ${rib} _3.tmp ${file}_ribAlign 2>&1| tee ${file}_rib.log   
bowtie -p3 -k2 -m1 -q --sam --un ${file%%.*}_ribUnalign ${rib} ${file%%.*}_SpikeInribUnalign ${file%%.*}_ribAlign 2>&1| tee ${file%%.*}_rib.log      
#cp -avr ${file}_ribAlign /home/macproadmin/users/Greg/Ribosom_alignFiles/
rm ${file%%.*}_ribAlign

#+++++++++++++++++++++++++++++++++++++++++++++++++++++HASHTAG IF NO rDNA+++++++++++++++++++++++++++++++++++++++++++++++++++++++++


for gen in ${genome}
do

echo ${file} aligning to genome...

#Generate 3 files: Unalign (not aligned to genome), Align (aligned to genome), _gen.log (bowtie stats)
#SAM format data is output from aligners that read FASTQ files and assign the sequences to a position with respect to a known reference genome.
bowtie -p3 -v2 -m1 -q --sam --un ${file%%.*}_Unalign ${gen} ${file%%.*}_ribUnalign ${file%%.*}_Align.sam 2>&1| tee ${file%%.*}_gen.log     


rm ${file%%.*}_ribUnalign
rm ${file%%.*}_Unalign

echo ${file%%.*} converting BAM to BED...
samtools view -bS ${file%%.*}_Align.sam > ${file%%.*}.bam
rm *.tmp
echo ${file%%.*} converting BAM to BED...
bamToBed -i ${file%%.*}.bam > ${file%%.*}.bed

echo indexing bam
samtools sort ${file%%.*}.bam ${file%%.*}_sorted
samtools index ${file%%.*}_sorted.bam ${file%%.*}_sorted.bai 
rm ${file%%.*}.bam


for chr in ${Combinedchrominfo}
do

#Converts BED to Sorted Bed using command line arguments
#sorts by chromosome
#Bedgraph file shows how many reads are present in each base position/
echo Sorting BED...
sort -k 1,1 -k2,2n ${file%%.*}.bed > ${file%%.*}_sorted.bed
rm ${file%%.*}.bed

# for splitting the sorted.bed file based on chromosomes for each organism

SC_chro=$(awk '{print $1}' ${SCchrominfo})
SP_chro=$(awk '{print $1}' ${SPchrominfo})


### IMPORTANT Note:  The steps below write repeatedly to the separated bed files.  Therefor, if you rerun this script without deleting or moving these files, they will be added to instead of overwritten. #######

for SCchr in ${SC_chro}
do
awk '$1=="'${SCchr}'" {print $0}' ${file%%.*}_sorted.bed >> ${file%%.*}_cerevisiae.bed 
sort -k 1,1 -k2,2n ${file%%.*}_cerevisiae.bed > ${file%%.*}_cerevisiae_sorted.bed
done
for SPchr in ${SP_chro}
do
awk '$1=="'${SPchr}'" {print $0}' ${file%%.*}_sorted.bed >> ${file%%.*}_pombe.bed
sort -k 1,1 -k2,2n ${file%%.*}_pombe.bed > ${file%%.*}_pombe_sorted.bed
done 

echo Number of reads mapping to S. cerevisiae genome: 
wc -l ${file%%.*}_cerevisiae_sorted.bed
echo Number of reads mapping to S. pombe genome: 
wc -l ${file%%.*}_pombe_sorted.bed
# remove unsorted bed files 
rm *e.bed

#Generate non-nomalized Bedgraphs 
#bedgraph - chr chrstart chrend

echo Generating non-normalized Bedgraphs pombe...
awk '$6 == "+"' ${file%%.*}_pombe_sorted.bed | genomeCoverageBed -i stdin -3 -bg -g ${SPchrominfo} > ${file%%.*}_pombe_plus.bedgraph 
awk '$6 == "-"' ${file%%.*}_pombe_sorted.bed | genomeCoverageBed -i stdin -3 -bg -g ${SPchrominfo} > ${file%%.*}_pombe_m.bedgraph
awk '{$4=$4*-1; print}' ${file%%.*}_pombe_m.bedgraph > ${file%%.*}_pombe_minus.bedgraph

#Generate normalized Bedgraphs for both organisms

echo Generate normalized Bedgraphs for pombe only
#set a normalization factor based on reads mapping to spike-In Genome reads per hundred thousand spike-in
nf=`wc ${file%%.*}_cerevisiae_sorted.bed | awk '{print $1}'` 
#nf=$(expr $spikeCounts/100000)
echo spikeCounts = $nf
#echo normalization factor is $nf
echo $nf | awk -v nf=$nf '{printf "%s\t%s\t%s\t%s\n", $1, $2, $3, ($4*100000/nf)}' ${file%%.*}_pombe_plus.bedgraph  > ${file%%.*}_pombe_SpikeNormed_plus.bedgraph
echo $nf | awk -v nf=$nf '{printf "%s\t%s\t%s\t%s\n", $1, $2, $3, ($4*100000/nf)}' ${file%%.*}_pombe_minus.bedgraph  > ${file%%.*}_pombe_SpikeNormed_minus.bedgraph

#Make BigWig
#Faster on browser, Used to generate further analysis
echo Making BigWigs...
${bgToBigWig} ${file%%.*}_pombe_plus.bedgraph ${SPchrominfo} ${file%%.*}_pombe_plus.bw
${bgToBigWig} ${file%%.*}_pombe_minus.bedgraph ${SPchrominfo} ${file%%.*}_pombe_minus.bw
${bgToBigWig} ${file%%.*}_pombe_SpikeNormed_plus.bedgraph ${SPchrominfo} ${file%%.*}_pombe_SpikeNormed_plus.bw
${bgToBigWig} ${file%%.*}_pombe_SpikeNormed_minus.bedgraph ${SPchrominfo} ${file%%.*}_pombe_SpikeNormed_minus.bw

#rm *.tmp
rm ${file%%.*}_pombe_m.bedgraph
rm *.sam
rm *Unalign*
gzip *.bed
gzip ${file%%.*}.fastq


echo 'RibAlign\tRepeatRib\tRepeatGen\tAlign\tUnalign\tTotal' > ${file}_stat.log 
echo `awk 'FNR == 2 {print $9}' ${file%%.*}_rib.log` "\t" `awk 'FNR == 4 {print $9}' ${file%%.*}_rib.log` "\t" `awk 'FNR == 4 {print $9}' ${file%%.*}_gen.log` "\t" `awk 'FNR == 2 {print $9}' ${file%%.*}_gen.log` "\t" `awk 'FNR == 3 {print $7}' ${file%%.*}_gen.log` "\t" `awk 'FNR == 1 {print $4}' ${file%%.*}_rib.log`>> ${file%%.*}_stat.log


#hashtag or remove one "done" if no rDNA
done
done
done
done
done


## check for existence of subdirectories, make if non-existant 
if [ ! -d "./bed" ]; then 
mkdir "./bed"
fi
if [ ! -d "./bedgraph" ]; then 
mkdir "./bedgraph"
fi
if [ ! -d "./logs" ]; then 
mkdir "./logs"
fi
if [ ! -d "./bw" ]; then 
mkdir "./bw"
fi
if [ ! -d "./bedgraph/bedgraph_FullRead" ]; then 
mkdir "./bedgraph/bedgraph_FullRead"
fi
if [ ! -d "./bw" ]; then 
mkdir "./bw"
fi
if [ ! -d "./bam" ]; then 
mkdir "./bam"
fi

## Move aligned files to corresponding directories
mv *.bed.gz ./bed 
mv *.bw ./bw
mv *.bedgraph ./bedgraph 
mv *.log ./logs
mv *.txt ./logs 
mv *.bam* ./bam
mv *.bai* ./bam
