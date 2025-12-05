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
4. cd to your parent folder and run ```source process.txt ``` 

## Features

- **Barcode Detection**: Uses LAST aligner to identify DNA barcodes in nanopore reads
- **Nucleosome Extraction**: Identifies and extracts nucleosome-sized fragments (13-1000 bp)
- **Quality Filtering**: Removes low-quality alignments (MAPQ < 40) and hairpin artifacts
- **Multiple Output Formats**: Generates BED, Juicer, and Cooler format files
- **Chromosome Conversion**: Converts chromosome names to Roman numerals for yeast genomes

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
3. **Chromosome sizes**: `sacCer2.chrsize3.txt`

# Pipeline Workflow

### **Part 1: common to all analysis, produce the `read_info.txt` file**

<ol start="1">
<li>
Index DNA barcodes:
<pre><code>lastdb Oligo_index/DNABarcode Oligo_index/DNABarcode.txt
</code></pre>
</li>
<li>
Concatenate and decompress FASTQ files:
<pre><code>cat fastq/*.fastq.gz > fastq/allreads.fastq.gz
gunzip fastq/allreads.fastq.gz
</code></pre>
</li>
<li>
Identify barcodes in reads:
<pre><code>cat fastq/allreads.fastq | parallel --pipe --recstart '\n@' -N10000 lastal -Q 1 Oligo_index/DNABarcode > out
perl SCRIPTS/maf_bcsplit_IW.pl out | grep -F + > maf_out_plus.txt
</code></pre>
</li>
<li>
Extract nucleosome positions (13–1000 bp):
<pre><code>sort -t',' -k1,1 -k4,4n maf_out_plus.txt | awk -F, -v OFS='\t' '{print $1,$4,$5,$7}' > maf_linkers.txt
mawk 'id != $1 {print id,pv2,pv2,pv2} {id=$1; pv2=$4}' maf_linkers.txt | sed 's/ /\t/g' >> maf_linkers1.txt
cat maf_linkers.txt maf_linkers1.txt > maf_linkers2.txt
sort -k1,1 -k2,2n maf_linkers2.txt > sortedMAF_linkers.txt
mawk 'id==$1 {print id,pv,$2,$2-pv} {id=$1; pv=$3}' sortedMAF_linkers.txt | mawk '$4 > 12 && $4 < 1001' | sed 's/ /\t/g' > maf_nucs.txt
</code></pre>
</li>
<li>
Create GTF and extract sequences:
<pre><code>mawk -v OFS='\t' '{print $1,"cad_C","nuc",$2,$3,$6=".",$7="+",$8=".",$9="."}' maf_nucs.txt > GTF_maf_nucs.txt
split -n 5 GTF_maf_nucs.txt
seqkit subseq --gtf xaa fastq/allreads.fastq > Maf_outa.fastq
seqkit subseq --gtf xab fastq/allreads.fastq > Maf_outb.fastq
seqkit subseq --gtf xac fastq/allreads.fastq > Maf_outc.fastq
seqkit subseq --gtf xad fastq/allreads.fastq > Maf_outd.fastq
seqkit subseq --gtf xae fastq/allreads.fastq > Maf_oute.fastq
cat Maf_out*.fastq > Maf_nucs.fastq
</code></pre>
</li>
<li>
Align to reference genome and convert to BED:
<pre><code>minimap2 -x map-ont -a --secondary=no -t 20 sacCer3.fa Maf_nucs.fastq > alignment.sam
samtools view -@16 -b alignment.sam > alignment.bam
bedtools bamtobed -i alignment.bam > alignment.bed
mawk '$5 > 40' alignment.bed | sed 's/_/\t/g' | sort -k4,4 -k1,1 -k5,5n > alignment_q40.bed
</code></pre>
</li>
<li>
Remove hairpin artifacts (±5 bp tolerance):
<pre><code>mawk -v OFS='\t' '((id>=($2-5) && id<=($2+5)) && rd==$4) || ((id2>=($3-5) && id2<=($2+5)) && rd==$4) {cond=0} (rd!=$4) {cond=1} {id=$2; id2=$3; rd=$4} { if (cond) print $1,$2,$3,$4,$5,$6,$7;}' alignment_q40.bed > alignment_dedup.bed
</code></pre>
</li>
<li>
Create read_info (homogenized format):
<pre><code>mawk '{gsub(/-/, "_", $4); print $1 "\t" $2 "\t" $3 "\t" $3-$2 "\t" $2+int(($3-$2)/2) "\t" $7 "\t" $4}' alignment_dedup.bed > read_info.txt
</code></pre>
</li>

This produces a tab separated text file 'read_info.txt' file with the following columns:
| Chromosome | Fragment start point | Fragment end point | Fragment size | Fragment midpoint | Strand| read_ID|
|------|------|------|------|------|------|-------------|
|chrXIII|679945|680261|316|680103|-|000000f8_4679_4275_a1c2_9d90bdd5a2d8 |
|chrXIII|538987|539124|137|539055|+|000000f8_4679_4275_a1c2_9d90bdd5a2d8 |
|chrXIII|549367|549479|112|549423|-|000000f8_4679_4275_a1c2_9d90bdd5a2d8 |
|chrXVI|412657|413085|428|412871|+|000000f8_4679_4275_a1c2_9d90bdd5a2d8 |
|chrXVI|413644|413775|131|413709|-|000000f8_4679_4275_a1c2_9d90bdd5a2d8 |

</ol>

### **Part 2: Produce high resolution map of interactions**

<ol start="9">
<li> Navigate to the `high_resolution_map` folder
<pre><code>cd high_resolution_map
source process.txt</code></pre>
</li>
<li>
Generate high-resolution pairs for Juicer:
<pre><code>python Readinfobed--pairs--CADC.py
sed "s/'//g" read_info.txt.juicer | mawk '{$1=$1; print}' > read_info.txt.sed.juicer
mawk '$2 > $6 {print $5,$6,$7,$8,$1,$2,$3,$4,$9} $2<=$6 {print}' read_info.txt.sed.juicer | parsort -k2,2d -k6,6d > read_info.s.juicer
</code></pre>
</li>
<li>
Generate contact matrices (.hic and .mcool):

**Replace here the path to your installation of** 
`juicer_tools.jar`

<pre><code>java -Xmx16g -jar /path/to/juicer_tools.jar pre read_info.s.juicer output.hic sacCer3 -r 10,25,33,50,75,100,150,200,300,400,500,1000,2500,5000,7500,10000,15000,20000,30000,50000,100000,250000,500000,1000000
mawk -F' ' '{print $2 "\t" $3 "\t" $6 "\t" $7 "\t"}' read_info.txt.sed.juicer > pairs
bedtools makewindows -g sacCer3.chrsize.txt -w 5 > bins.5.txt
cooler cload pairs -c1 1 -p1 2 -c2 3 -p2 4 bins.5.txt pairs extr.cool
cooler zoomify -r 5N extr.cool
</code></pre>
</li>
<li>
Sort read_info by genomic coordinates:
<pre><code>parsort -k1,1 -k2,2n read_info.txt > read_info.r.s.bed
</code></pre>
</li>
</ol>

### **Part 3: Produce map of sister nucleosome ligation**

<ol start="13">
<li>
Generate read_info-based Juicer pairs and normalize IDs:
<pre><code>cd ../sister_chromatid_overlap
source process.txt
</code></pre>
</li>

<li>
Generate read_info-based Juicer pairs and normalize IDs:
<pre><code>python Readinfobed--pairs--CADC_test.py
sed "s/'//g" read_info.txt_test.juicer | mawk '{$1=$1; print}' > read_info.txt.sed.juicer
mawk '$2 > $6 {print $5,$6,$7,$8,$1,$2,$3,$4,$9} $2<=$6 {print}' read_info.txt.sed.juicer | parsort -k2,2d -k6,6d > read_info.s.juicer
</code></pre>
</li>
<li>
Build pairwise intervals (around midpoints) and clean:
<pre><code>mawk '{print $2 "\t" $3-int($1/2) "\t" $3+int($1/2) "\t" $6 "\t" $7-int($5/2) "\t" $7+int($5/2) "\t" $9 }' read_info.s.juicer > pairwise.bed
rg -v -e '2micron' -e 'chrM' pairwise.bed > pairwise.cln.bed
bedtools overlap -i pairwise.cln.bed -cols 2,3,5,6 > pairwise.olp.bed
</code></pre>
</li>
<li>
Select overlaps on same contig with overlap length > 3:
<pre><code>mawk '$1 == $4 && $8 > 3' pairwise.olp.bed > allolp.bed
</code></pre>
</li>
<li>
Add length and ratio columns to quantify overlaps:
<pre><code>mawk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $3-$2 "\t" $6-$5 "\t" int((($3-$2)+($6-$5))/2) "\t"  int((($3-$2)+($6-$5))/2)/$8 }' allolp.bed > olp.pct.bed
</code></pre>
</li>
<li>
Filter overlaps by ratio, minimum overlap length, and region size ranges:
<pre><code>mawk '(($12 <= 0.98 || $12 >= 1.02) && $8 > 15 && $9 >= 100 && $9 <= 1000 && $10 >= 100 && $10 <= 1000)' olp.pct.bed > olp.pct.filt.bed
</code></pre>
</li>
<li>
Extract actual overlap coordinates (max of starts, min of ends):
<pre><code>mawk '{start = ($2 > $5 ? $2 : $5); end = ($3 < $6 ? $3 : $6); print $1 "\t" start "\t" end}' olp.pct.filt.bed > real-olp4.bed
</code></pre>
</li>
<li>
Sort, filter out positions ≤1, and convert to bedGraph:
<pre><code>parsort -k1,1 real-olp4.bed > real-olp4.s.bed
mawk '($2>1)' real-olp4.s.bed > real-olp4.sf.bed
bedtools genomecov -i real-olp4.sf.bed -g ../sacCer3.chrsize.txt -bg > real-olp4.bedgraph
</code></pre>
</li>
<li>
Generate BigWig for raw overlap coverage:
<pre><code>bedGraphToBigWig real-olp4.bedgraph ../sacCer3.chrsize.txt real-olp4.bigwig
</code></pre>
</li>
<li>
Compute CPM-normalized bedGraph and BigWig: 

**Replace here the number by your number of reads as detected by** 
`wc -l read_info.txt`
<pre><code>wc -l read_info.txt
mawk '{print $1 "\t" $2 "\t" $3 "\t" $4/195}' real-olp4.bedgraph > real-olp4.cpm.bedgraph
bedGraphToBigWig real-olp4.cpm.bedgraph ../sacCer3.chrsize.txt real-olp4.cpm.bigwig
</code></pre>
</li>
<li>
Prepare per-read records for the filtered overlaps:
<pre><code>mawk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7}' olp.pct.filt.bed > olp.pct.filt1.bed
mawk 'NR==FNR {ids[$7]; next} ($6 in ids)' olp.pct.filt.bed read_info.txt > olp.read_info.txt
python Readinfobed--pairs--CADC_OLP.py
sed "s/'//g" olp.read_info.txt_test.juicer | mawk '{$1=$1; print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8}' > olp.read_info.txt.sed.juicer
</code></pre>
</li>
<li>
Order Juicer pairs and generate `.hic`:
<pre><code>mawk '$2 > $6 {print $5,$6,$7,$8,$1,$2,$3,$4,$9} $2<=$6 {print}' olp.read_info.txt.sed.juicer | parsort -k2,2d -k6,6d > olp.read_info.s.juicer
java -Xmx16g -jar /path/to/juicer_tools.jar pre olp.read_info.s.juicer olp.hic sacCer3 -r 10,25,33,50,75,100,150,200,300,400,500,1000,2500,5000,7500,10000,15000,20000,30000,50000,100000,250000,500000,1000000
</code></pre>
</li>
<li>
Create Cooler file from pairs and zoomify:
<pre><code>mawk -F' ' '{print $2 "\t" $3 "\t" $6 "\t" $7 "\t"}'  olp.read_info.txt.sed.juicer > olp.pairs
bedtools makewindows -g ../sacCer3.chrsize.txt -w 5 > bins.5.txt
cooler cload pairs -c1 1 -p1 2 -c2 3 -p2 4 bins.5.txt olp.pairs olp.extr.cool
cooler zoomify -r 5N olp.extr.cool
</code></pre>
</li>
<li>
Counts: quantify long reads and pair totals:
<pre><code>mawk -F'\t' '{count[$6]++} END {for (val in count) print val, count[val]}' olp.read_info.txt > olp.longread.txt
mawk -F'\t' '{count[$7]++} END {for (val in count) print val, count[val]}' ../read_info.txt > all.longread.txt
wc -l olp.longread.txt all.longread.txt olp.pairs ../pairs > counts.olp.longread.txt
</code></pre>
</li>
</ol>

## Output Files

| File | Description |
|------|-------------|
| read_info.txt | Processed read information |
| *.hic | Juicer Hi-C format file |
| *.mcool | Cooler format contact matrix |
|*.olp.cpm.bigwig | BigWig track of sister chromatid 
|*.olp.hic | Olp + 1 interaction for Juicebox
|*.olp.mcool |Olp + 1 interaction for HiGlass
## Directory Structure

```bash
├── fastq/                  # Input FASTQ files
├── Oligo_index/           # DNA barcode sequences
│   └── DNABarcode.txt
├── SCRIPTS/               # Helper scripts
│   ├── maf_bcsplit_IW.pl
│   └── Readinfobed--pairs--CADC.py
├── high_resolution_map       
│   ├── *.hic               # CAD-C map for Juicebox
│   └── *.mcool             # CAD-C map for HiGlass
├── sister_crhomatid_overlap       
│   ├── *.olp.cpm.bigwig        # BigWig track of sister chromatid nucleosome overlap (Olp)
│   ├── *.olp.hic               # Olp + 1 interaction for Juicebox
│   └── *.olp.mcool             # Olp + 1 interaction for HiGlass
└── README.md
```



