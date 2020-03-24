#PBS -S /bin/bash
#PBS -N Dap-SeqPipeline
#PBS -q batch
#PBS -l nodes=1:ppn=24
#PBS -l walltime=120:00:00
#PBS -l mem=40gb

## should be prepared..
# 1. Data // 0.Reference 

## What should I change for that
#These!!..
Data_Path=(/scratch/sb14489/2-6.DAPseq_soybean_v4_withGST/1.Data)
Path=(/scratch/sb14489/2-6.DAPseq_soybean_v4_withGST)
RawDataFormat=.fastq.gz
ControlSample=GST_Soybean  #if you don't have any control then put "NA"
Reference=Gmax_508_v4.0_OnlyChr #Don't add .fa or path

########



##### Start! ###########
List=`find "$Data_Path" -name "*""$RawDataFormat" | sed 's|.*/||'`

cd $Path

echo "This is input fastq!: "$List

#function FASTAQC_forRawData(){
#}

function Trimmomatic(){
module load Trimmomatic/0.36-Java-1.8.0_144

mkdir -p  2.Trimmomatic
for i in $List;
do 
item="${i/$RawDataFormat/}"
#echo $item
#<<"END"
echo java -jar /usr/local/apps/eb/Trimmomatic/0.36-Java-1.8.0_144/trimmomatic-0.36.jar \
       SE -threads 24  -phred33 "$Data_Path"/"$i" ./2.Trimmomatic/"$item"_trimmed.fastq.gz \
       ILLUMINACLIP:/scratch/sb14489/Program/trimmomatic-0.36/adapters/TruSeq3-SE.fa:2:30:10\
       LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50  2> ./2.Trimmomatic/"$item".log

java -jar /usr/local/apps/eb/Trimmomatic/0.36-Java-1.8.0_144/trimmomatic-0.36.jar \
       SE -threads 24  -phred33 "$Data_Path"/"$i" ./2.Trimmomatic/"$item"_trimmed.fastq.gz \
       ILLUMINACLIP:/scratch/sb14489/Program/trimmomatic-0.36/adapters/TruSeq3-SE.fa:2:30:10\
       LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50  2> ./2.Trimmomatic/"$item".log

#END
done
}

#function FASTAQC_forTrimmedData(){
#}

########################################################################################
function BuildFASTA(){

module load Bowtie2/2.3.4.1-foss-2016b
#bowtie2-build --threads 24 ./0.Reference/"$Reference".fa  ./0.Reference/"$Reference"

module load SAMtools/1.9-foss-2016b
samtools faidx ./0.Reference/"$Reference".fa
cut -f1,2 ./0.Reference/"$Reference".fa.fai > ./0.Reference/"$Reference".chrom.sizes
ml pyfaidx/0.5.5.1-foss-2016b-Python-2.7.14

mkdir -p  ./0.Reference/"$Reference"_chromosomes
cp ./0.Reference/"$Reference".fa ./0.Reference/"$Reference"_chromosomes/ 
cd ./0.Reference/"$Reference"_chromosomes/
faidx -x "$Reference".fa 
}
########################################################################################
function Bowtie2(){
cd $Path
mkdir -p 3.Bowtie2
module load Bowtie2/2.3.5.1-GCC-8.2.0-2.31.1

for i in $List;
do
item="${i/$RawDataFormat/}"
Outfile=("$Path""/3.Bowtie2/""$item"".sam")
Logfile=("$Path""/3.Bowtie2/""$item"".log")

bowtie2 -x "$Path"/0.Reference/"$Reference" -U ./2.Trimmomatic/"$item"_trimmed.fastq.gz -S $Outfile -p 24 2> $Logfile

done

}
########################################################################################
function Sorting_bg_bw(){
module load SAMtools/1.10-GCC-8.2.0-2.31.1
module load BEDTools/2.28.0-foss-2018a
module load ucsc/359

cd "$Path"
mkdir -p 4.Sorted_Filtered_BamFile
for i in $List;
do
item="${i/$RawDataFormat/}"
samtools sort -O 'bam' -o ./4.Sorted_Filtered_BamFile/${item}_sorted.bam ./3.Bowtie2/${item}.sam
samtools view -q 30 ./4.Sorted_Filtered_BamFile/${item}_sorted.bam -o ./4.Sorted_Filtered_BamFile/${item}.sorted.filtered.bam
samtools index ./4.Sorted_Filtered_BamFile/${item}.sorted.filtered.bam
bedtools bamtobed -i ./4.Sorted_Filtered_BamFile/${item}.sorted.filtered.bam > ./4.Sorted_Filtered_BamFile/${item}.bed

bedtools genomecov -i ./4.Sorted_Filtered_BamFile/${item}.bed -split -bg -g ./0.Reference/"$Reference".chrom.sizes  > ./4.Sorted_Filtered_BamFile/${item}.bg
wigToBigWig ./4.Sorted_Filtered_BamFile/${item}.bg ./0.Reference/"$Reference".chrom.sizes  ./4.Sorted_Filtered_BamFile/${item}.bw

done
}

########################################################################################
function GEM(){
cd "$Path"
mkdir -p 5.GEM

if [ "$ControlSample" = "NA" ];
then
for i in $List;
do
item="${i/$RawDataFormat/}"
java -Xmx20G -jar /scratch/sb14489/Program/gem/gem.jar --d /scratch/sb14489/Program/gem/Read_Distribution_default.txt --g ./0.Reference/${Reference}.chrom.sizes --genome ./0.Reference/${Reference}_chromosomes --f BED  --expt ./4.Sorted_Filtered_BamFile/${item}.bed --out ./5.GEM/${item} --t 24 --k_min 6 --k_max 20 --outNP --sl --q 5
done
else
NewList=( "${List[@]/$ControlSample}" )
for i in $NewList;
do
item="${i/$RawDataFormat/}"
java -Xmx20G -jar /scratch/sb14489/Program/gem/gem.jar --d /scratch/sb14489/Program/gem/Read_Distribution_default.txt --g ./0.Reference/${Reference}.chrom.sizes --genome ./0.Reference/${Reference}_chromosomes --f BED  --expt ./4.Sorted_Filtered_BamFile/${item}.bed --ctrl ./4.Sorted_Filtered_BamFile/${ControlSample}.bed --out ./5.GEM/${item} --t 24 --k_min 6 --k_max 20 --outNP --sl --q 5

done
fi
}
#######################################################################################
function Counting_PeakNumber(){
python /scratch/sb14489/Program/AfterGEM_CountingPeakNumb_PrepareChIPQC.py
}

#Trimmomatic
#BuildFASTA
#Bowtie2
#Sorting_bg_bw
#GEM
Counting_PeakNumber
