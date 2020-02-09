#https://www.ncbi.nlm.nih.gov/sra/docs/sradownload/
#using prefetch https://github.com/ncbi/sra-tools/wiki/Download-On-Demand

module load SRA-Toolkit/2.9.1-centos_linux64


List=(SRX3786361 SRX3786362 SRX3786363 SRX3786364 SRX3786365 SRX3786366 SRX3786367 SRX3786368 SRX3786369 SRX3786370 SRX3786371 SRX3786372 SRX3786373 SRX3786374)
for item in ${List[@]};
do
echo prefetch "$item" -O  /scratch/sb14489/2.DAPseq/1.DownloadData/
prefetch "$item" -O  /scratch/sb14489/2.DAPseq/1.DownloadData/
done

