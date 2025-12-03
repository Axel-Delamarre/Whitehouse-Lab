# CAD-C Data Processing Pipeline

A bioinformatics pipeline for processing CAD-C data from Oxford Nanopore sequencing to identify and map nucleosome interactions.

## Overview

This pipeline processes nanopore sequencing reads to:
- Identify DNA barcodes/linkers in long reads
- Extract individual nucleosome sequences while retaining read connectivity information
- Map nucleosomes to reference genome
- Generate contact matrices for chromatin interaction analysis
- Produce Hi-C compatible output files

## To run this pipeline:

1. Install the required tools  
2. Copy this whole directory on your workspace
3. Place your ```.fastq``` in the fastq folder
4. cd to your parent folder and run ```source source.txt ``` 

## Features

- **Barcode Detection**: Uses LAST aligner to identify DNA barcodes in nanopore reads
- **Nucleosome Extraction**: Identifies and extracts nucleosome-sized fragments (13-1000 bp)
- **Quality Filtering**: Removes low-quality alignments (MAPQ < 40) and hairpin artifacts
- **Multiple Output Formats**: Generates BED, Juicer, and Cooler format files

## Prerequisites
### Core tools
1. LAST Aligner : https://gitlab.com/mcfrith/last 
2. minimap2 : https://github.com/lh3/minimap2 
3. samtools : http://www.htslib.org/ 
4. bedtools : https://bedtools.readthedocs.io/ 
5. seqkit : https://bioinf.shenwei.me/seqkit/ 
6. GNU parallel : https://www.gnu.org/software/parallel/ 
7. mawk : https://invisible-island.net/mawk/ GNU awk : https://www.gnu.org/software/gawk/ 
8. Perl : https://www.perl.org/ 
9. cooler : https://cooler.readthedocs.io/ 

### Python packages
- Python 3.x
- Custom script: Readinfobed--pairs--CADC.py

### Java tools
- Juicer tools (juicer_tools.jar)

## Input Files Required

1. **Nanopore FASTQ files**: Place in `fastq/` directory (`.fastq.gz` format)
2. **DNA Barcode file**: `Oligo_index/DNABarcode.txt`
3. **Reference genome**: `sacCer3.fa`
4. **Chromosome sizes**: `sacCer2.chrsize3.txt`

## Pipeline Workflow

- **Step 1: Index DNA Barcodes**
```bash
lastdb Oligo_index/DNABarcode.txt Oligo_index/DNABarcode.txt
```
- **Step 2: Concatenate and Decompress FASTQ Files**
```bash
cat fastq/*.fastq.gz > fastq/allreads.fastq.gz
gunzip fastq/allreads.fastq.gz
```
- **Step 3: Identify Barcodes in Reads**
```bash
cat Fastq/*.fastq | parallel --pipe --recstart '\n@' -N10000 lastal -Q 1 Oligo_index/DNABarcode.txt > out
perl SCRIPTS/maf_bcsplit_IW.pl out | grep + > maf_out_plus.txt
```
- **Step 4: Extract Nucleosome Positions**

Sort and process linker positions
```bash
sort -t',' -k1,1 -k4,4n maf_out_plus.txt | awk -F, -v OFS='\t' '{print $1,$4,$5,$7}' > maf_linkers.txt
```
Insert extra line to capture read ends
```bash
mawk 'id != $1 {print id,pv2,pv2,pv2} {id=$1; pv2=$4}' maf_linkers.txt | sed 's/ /\t/g' >> maf_linkers1.txt
```
Combine and sort
```bash
cat maf_linkers.txt maf_linkers1.txt > maf_linkers2.txt
sort -k1,1 -k2,2n maf_linkers2.txt > sortedMAF_linkers.txt
```
Calculate insert sizes (13-1000 bp range)
```bash
mawk 'id==$1 {print id,pv,$2,$2-pv} {id=$1; pv=$3}' sortedMAF_linkers.txt | mawk '$4 > 12 && $4 < 1001' | sed 's/ /\t/g' > maf_nucs.txt
```
- **Step 5: Create GTF and Extract Sequences**

Create GTF file
```bash
mawk -v OFS='\t' '{print $1,"cad_C","nuc",$2,$3,$6=".",$7="+",$8=".",$9="."}' maf_nucs.txt > GTF_maf_nucs.txt
```
Split GTF file for parallel processing
```bash
split -n 5 GTF_maf_nucs.txt
```

Extract sequences using seqkit
```bash
seqkit subseq --gtf xaa fastq/allreads.fastq > Maf_outa.fastq
seqkit subseq --gtf xab fastq/allreads.fastq > Maf_outb.fastq
seqkit subseq --gtf xac fastq/allreads.fastq > Maf_outc.fastq
seqkit subseq --gtf xad fastq/allreads.fastq > Maf_outd.fastq
seqkit subseq --gtf xae fastq/allreads.fastq > Maf_oute.fastq
```
Concatenate results
```bash
cat Maf_out*.fastq > Maf_nucs.fastq
```
- **Step 6: Align to Reference Genome**

Map with minimap2
```bash
minimap2 -amap-ont --secondary=no --sam-hit-only -t20 sacCer3.fa Maf_nucs.fastq > alignment.sam
```

Convert to BED (for large files, use BAM method)
```bash
samtools view -@16 -b -S alignment.sam > file.bam
bamtobed -i file.bam > alignment.bed
```
Filter by quality (MAPQ > 40)
```bash
cat alignment.bed | awk '$5 >40' | sed 's/_/\t/g' | sort -k4,4 -k1,1 -k5,5n > alignment_q40.bed
```
- **Step 7: Remove Hairpin Artifacts**

Remove duplicate reads from hairpin structures (±5bp tolerance)
```bash
mawk -v OFS='\t' '((id>=($2-5) && id<=($2+5)) && rd==$4) || ((id2>=($3-5) && id2<=($2+5)) && rd==$4) {cond=0} (rd!=$4) {cond=1} {id=$2; id2=$3; rd=$4} { if (cond) print $1,$2,$3,$4,$5,$6,$7;}' alignment_q40.bed > alignment_dedup.bed
```
- **Step 8: Process Read Information**

Create read info file
```bash
awk '{gsub(/-/, "_", $4); print $1 "\t" $2 "\t" $3 "\t" $3-$2 "\t" $2+int(($3-$2)/2) "\t" $7 "\t" $4}' alignment_dedup.bed > read_info.txt
```
Generate pairs for Juicer
```bash
python Readinfobed--pairs--CADC.py
```
Clean up formatting
```bash
sed "s/'//g" read_info.txt.juicer | awk '{$1=$1; print}' > read_info.txt.sed.juicer
```
Sort for Juicer
```bash
awk '$2 > $6 {print $5,$6,$7,$8,$1,$2,$3,$4,$9} $2<=$6 {print}' read_info.txt.sed.juicer | parsort -k2,2d -k6,6d > read_info.s.juicer
```
- **Step 9: Generate Contact Matrices**
Juicer Hi-C Format
```bash
java -Xmx16g -jar /path/to/juicer_tools.jar pre read_info.s.juicer output.hic sacCer3 -r 10,25,33,50,75,100,150,200,300,400,500,1000,2500,5000,7500,10000,15000,20000,30000,50000,100000,250000,500000,1000000
```
Cooler Format
```bash
awk -F' ' '{print $2 "\t" $3 "\t" $6 "\t" $7 "\t"}' read_info.txt.sed.juicer > pairs
bedtools makewindows -g sacCer2.chrsize3.txt -w 5 > bins.5.txt
cooler cload pairs -c1 1 -p1 2 -c2 3 -p2 4 bins.5.txt pairs extr.cool
cooler zoomify -r 5N extr.cool
```
- **Step 10: Sort the read_info file**
```bash
parsort -k1,1 -k2,2n read_info.txt > read_info.r.s.bed
```


## Output Files

| File | Description |
|------|-------------|
| maf_nucs.txt | Nucleosome positions between linkers |
| read_ID_countsq.txt | Frequency of nucleosomes per read |
| GTF_maf_nucs.txt | GTF format file of nucleosome positions |
| Maf_nucs.fastq | Extracted nucleosome sequences |
| alignment_dedup.bed | Deduplicated alignment coordinates |
| read_info.txt | Processed read information |
| *.hic | Juicer Hi-C format file |
| extr.cool | Cooler format contact matrix |
| read_info.r.s.bed | Sorted BED with Roman numeral chromosomes |

## Parameters

### Key Thresholds
- **Insert size range**: 13-1000 bp (nucleosome-sized fragments)
- **Mapping quality**: MAPQ > 40
- **Hairpin detection tolerance**: ±5 bp
- **Alignment mode**: Oxford Nanopore (`-amap-ont`)
- **Threading**: 20 threads (adjust with `-t` flag)


## Directory Structure

```bash
├── fastq/                  # Input FASTQ files
├── Oligo_index/           # DNA barcode sequences
│   └── DNABarcode.txt
├── SCRIPTS/               # Helper scripts
│   ├── maf_bcsplit_IW.pl
│   └── Readinfobed--pairs--CADC.py
├── sacCer3.fa            # Reference genome
├── sacCer2.chrsize3.txt  # Chromosome sizes
└── README.md
```


