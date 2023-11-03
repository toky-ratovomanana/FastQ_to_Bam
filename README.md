# FastQ_to_Bam

To switch from fastQ (R1 & R2) to a BAM, you need to do the mapping with BWA and create the bam with Samtools.
To do this, you need the reference genome and and its index (Here: Homo_sapiens.GRCh38.ensembl_91.dna.primary_assembly.fa).

### Installation

Install BWA: https://installati.one/install-bwa-ubuntu-20-04/
```
sudo apt-get update
sudo apt-get -y install bwa
```

Install Samtools: https://www.htslib.org/download/
```
cd ..
wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2
tar -vxjf samtools-1.9.tar.bz2
cd samtools-1.9
make
```

### Example
FastQ_R1: CCARE_SA_075No_S17_R1_001.fastq.gz
FastQ_R2: CCARE_SA_075No_S17_R2_001.fastq.gz

Name of Bam output: CCARE_SA_075No_S17.bam
```
bwa mem -K 100000000 -v 3 -t 4 Homo_sapiens_ensembl38/Homo_sapiens.GRCh38.ensembl_91.dna.primary_assembly.fa <(zcat Fastq/CCARE_SA_075No_S17_R1_001.fastq.gz) <(zcat Fastq/CCARE_SA_075No_S17_R2_001.fastq.gz) | samtools view -b -o CCARE_SA_075No_S17.bam
```
