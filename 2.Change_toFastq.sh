#PBS -S /bin/bash
#PBS -N Convert_fastq
#PBS -q batch
#PBS -l nodes=1:ppn=4
#PBS -l walltime=120:00:00
#PBS -l mem=20gb

module load SRA-Toolkit/2.9.1-centos_linux64

cd /scratch/sb14489/2.DAPseq/1.DownloadData

List=(SRX3786361 SRX3786362 SRX3786363 SRX3786364 SRX3786365 SRX3786366 SRX3786367 SRX3786368 SRX3786369 SRX3786370 SRX3786371 SRX3786372 SRX3786373 SRX3786374)
for item in ${List[@]};
do
echo fastq-dump "$item" 2> 2.ChangeName.log
fastq-dump "$item" 2> 2.ChangeName.log
done
