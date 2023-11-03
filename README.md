# FastQ_to_Bam


### Installation

Install BWA: https://installati.one/install-bwa-ubuntu-20-04/
```
sudo apt-get update
sudo apt-get -y install bwa
```

### Example
FastQ_R1: CCARE_SA_075No_S17_R1_001.fastq.gz
FastQ_R2: CCARE_SA_075No_S17_R2_001.fastq.gz

Name of Bam output: CCARE_SA_075No_S17.bam
```
bwa mem -K 100000000 -v 3 -t 4 Homo_sapiens_ensembl38/Homo_sapiens.GRCh38.ensembl_91.dna.primary_assembly.fa <(zcat Fastq/CCARE_SA_075No_S17_R1_001.fastq.gz) <(zcat Fastq/CCARE_SA_075No_S17_R2_001.fastq.gz) | samtools view -b -o CCARE_SA_075No_S17.bam
```
