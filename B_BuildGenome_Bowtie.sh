#PBS -S /bin/bash
#PBS -N Convert_fastq
#PBS -q batch
#PBS -l nodes=1:ppn=16
#PBS -l walltime=120:00:00
#PBS -l mem=10gb


module load Bowtie2/2.3.5.1-GCC-8.2.0-2.31.1

cd /scratch/sb14489/2.DAPseq/0.Reference/

bowtie2-build --threads 16 Zea_mays.AGPv4.dna.toplevel_OnlyChr.fa Zea_mays.AGPv4.dna.toplevel_OnlyChr
