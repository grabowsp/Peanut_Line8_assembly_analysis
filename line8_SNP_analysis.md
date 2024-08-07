# Scripts used for SNP analysis of Line8

## Requirements
* bwa-mem2
* samtools
* gzip
* picard
* varscan
* bcftools
* SnpEff

## Steps
1. Index Tifrunner assembly with bwa and samtools
2. Align Line8 and parent Illumina reads to Tifrunner with bwa-mem2
3. Remove duplicate reads with picard
4. Generate gVCF with varscan
5. Estimate effects of SNPs with SnpEff

## Input files
Includes file names used in scripts below
* Tifrunner reference files
  * Tifrunner fasta: Tifrunner.fa
  * Tifrunner genes gff3 file: Tifrunner.gff3
  * Tifrunner CDS fasta: Tifrunner.cds.fa
  * Tifrunner protein fasta: Tifrunner.faa
* Fastq files of Illumina reads
  * Line8: Line8.fastq
  * C76-16: C76.fastq
  * Georgia Green: GG.fastq

## Index Tifrunner assembly
```
bwa-mem2 index Tifrunner.fa
samtools faidx Tifrunner.fa
```

## Align Illumina reads to Tifrunner
```
bash

IN_FILE_PREFIX=Line8 	# repeat following steps using C76 and GG as prefix
FASTQ_FILE=$IN_FILE_PREFIX'.fastq'
TMP_DIR=/temporary/directory
OUT_BAM=$IN_FILE_PREFIX'.bam'

bwa-mem2 mem -M -t8 -p Tifrunner.fa $FASTQ_FILE | \
samtools view -u -f 2 -F 4 - | \
samtools sort -@ 4 -T $TMP_DIR -o $OUT_BAM -
```

## Remove duplicate reads
```
bash

IN_FILE_PREFIX=Line8    # repeat following steps using C76 and GG as prefix
IN_BAM=$IN_FILE_PREFIX'.bam'
OUT_BAM=$IN_FILE_PREFIX'.dedup.bam'
OUT_STATS=$IN_FILE_PREFIX'.stats'
TMP_DIR=/temporary/directory

picard MarkDuplicates I=$IN_BAM O=$OUT_BAM M=$OUT_STATS VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true SORTING_COLLECTION_SIZE_RATIO=0.05 REMOVE_DUPLICATES=true TMP_DIR=$TMP_DIR MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 MAX_RECORDS_IN_RAM=2750000 READ_NAME_REGEX='[a-zA-Z0-9]+[:_][0-9]+[:_][a-zA-Z0-9]+[:_][0-9]+[:_]([0-9]+)[:_]([0-9]+)[:_]([0-9]+)' CREATE_INDEX=true

```

## Make gVCFs
### Header text for final gVCFs
* Adjust the last field of the bottom line to the correct sample
* Adjust 'fileDate' accordingly
* save as `header.txt`
```
##fileformat=VCFv4.2
##fileDate=DayMonthYear
##reference=Tifrunner.fa
##source=VarScan2cns
##FILTER=<ID=str10,Description="Less than 10% or more than 90% of variant supporting reads on one strand">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=RD,Number=1,Type=Integer,Description="Depth of reference-supporting bases (reads1)">
##FORMAT=<ID=AD,Number=1,Type=Integer,Description="Depth of variant-supporting bases (reads2)">
##FORMAT=<ID=PVAL,Number=1,Type=String,Description="P-value from Fisher's Exact Test">
#CHROM     	POS 	ID     	REF  	ALT  	QUAL   FILTER  INFO	FORMAT      	Line8
```
### Make gVCF from deduplicated bam
```
bash

IN_FILE_PREFIX=Line8    # repeat following steps using C76 and GG as prefix
IN_BAM=$IN_FILE_PREFIX'.dedup.bam'
TMP_VCF=$IN_FILE_PREFIX'.tmp.vcf.gz'
FULL_VCF=$IN_FILE_PREFIX'.vcf.gz'
DEPTH8_VCF=$IN_FILE_PREFIX'.d8.vcf.gz'

samtools mpileup -d 500 -Q 20 -f Tifrunner.fa $IN_BAM | \
varscan mpileup2cns --min-coverage 1 --min-reads2 0 --output-vcf | \
grep -v "#" | \
awk 'FS=OFS="\t" {split($10,s,":"); if(s[1] != "./.") print $1,$2,$3,$4,$5,$6,$7,".","GT:RD:AD", s[1]":"s[5]":"s[6]}' | \
bgzip -c > $TMP_VCF;

bgzip -c header.txt > header.txt.gz
cat header.txt.gz $TMP_VCF > $FULL_VCF
bcftools index $FULL_VCF

# Filter by sites with 8 or more reads
gzip -dck $FULL_VCF | \
awk '/^#/ {print} OFS="\t" {split($10,a,":"); split(a[2],b,","); $10 = a[1]; $9="GT"; if(b[1] + b[2] >= 8) print $0}'  | \
bgzip -c > DEPTH8_VCF

```

## Run SnpEff on Line8 VCF
### Adjust config file
```
# copy the standard config file to Tifrunner.config, then add following 
#   text to bottom of file

# Tifrunner peanut genome
arahy.Tifrunner.genome : Tifunner

<param name="genomeVersion" type="select" label="Genome">
    <option value="arahy.Tifrunner">Tifrunner (arahy.Tifrunner)<option>
<param>
```
### Prepare directory to make SnpEff database
```
# make subdirectory in directory with Tifrunner.config file
mkdir arahy.Tifrunner
cd arahy.Tifrunner

cp Tifrunner.fa ./sequences.fa
cp Tifrunner.gff3 ./genes.gff
cp Tifrunner.cds.fa ./cds.fa
cp Tifrunner.faa ./protein.fa
```
### Make SnpEff database
* run from directory with Tifrunner.config file
```
snpEff build -gff3 -v arahy.Tifrunner -config ./Tifrunner.config -dataDir .
```
### Get SnpEff annotations for VCF
```
bash

# copy Line8.d8.vcf.gz to directory with Tifrunner.config file

gunzip Line8.d8.vcf.gz

snpEff arahy.Tifrunner Line8.d8.vcf.gz -config ./Tifrunner.config -datadir . > \
Line8.d8.annotated.vcf

```
