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
> Tester sur mon PC avec Version: 1.10 (using htslib 1.10.2-3ubuntu0.1)
```
cd ..
wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2
tar -vxjf samtools-1.9.tar.bz2
cd samtools-1.9
make
```

### Example (Run Mapping)
- FastQ_R1: CCARE_SA_075No_S17_R1_001.fastq.gz
- FastQ_R2: CCARE_SA_075No_S17_R2_001.fastq.gz

Name of Bam output: CCARE_SA_075No_S17.bam
```
bwa mem -K 100000000 -v 3 -t 4 Homo_sapiens_ensembl38/Homo_sapiens.GRCh38.ensembl_91.dna.primary_assembly.fa <(zcat Fastq_batch_1/CCARE_SA_075No_S17_R1_001.fastq.gz) <(zcat Fastq_batch_1/CCARE_SA_075No_S17_R2_001.fastq.gz) | samtools view -b -o CCARE_SA_075No_S17.bam
```

### To run multiple files you can use nextflow script
First Install Nextflow: https://www.nextflow.io/docs/latest/getstarted.html
```
curl get.nextflow.io | bash
sudo mv nextflow /usr/local/bin
nextflow self-update
```
#To run the nextflow script (wes_map_sort_option_trim.nf)
* --inputDir containg the fastq filles (R1 & R2)
* --outDir directory created for output
* --refGenome Directory with the all the reference genome file (here 44Go)

```
NXF_VER=22.04.5 nextflow run wes_map_sort_option_trim.nf --inputDir ~/Documents/07_MSInsight/05_Work/11_Hokla/02_Dev-Frontend_Backend/fastQ_to_Bam/Fastq_batch_1 --outDir ~/Documents/07_MSInsight/05_Work/11_Hokla/02_Dev-Frontend_Backend/fastQ_to_Bam/Mapping_ensembl_38 --refGenome ~/Documents/07_MSInsight/05_Work/11_Hokla/02_Dev-Frontend_Backend/fastQ_to_Bam/Homo_sapiens_ensembl38/Homo_sapiens.GRCh38.ensembl_91.dna.primary_assembly.fa --trim F
```
