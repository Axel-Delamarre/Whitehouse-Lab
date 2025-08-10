## PCP Data Processing Workflow

### Tool Requirements

- **UMI-tools**: [Documentation](https://umi-tools.readthedocs.io/en/latest/)
- **seqkit**: [Documentation](https://bioinf.shenwei.me/seqkit/)
- **trim_galore**: [GitHub](https://github.com/FelixKrueger/TrimGalore)
- **bowtie2**: [bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
- **samtools**: [samtools](https://www.htslib.org/)
- **python**: [python](https://www.python.org/downloads/)
- **juicer_tools**: [GitHub](https://github.com/aidenlab/JuicerTools)
- **cooler**: [GitHub](https://github.com/open2c/cooler)
- **bedtools**: [bedtools](https://github.com/arq5x/bedtools2)

**Visualization:**

- `.hic` files: [Juicebox](https://github.com/aidenlab/Juicebox)
- `.mcool` files: [HiGlass](https://higlass.io/) or [pyGenomeTracks](https://pygenometracks.readthedocs.io/en/latest/)

---

### **Preparation**

Before proceeding, ensure all required programs are installed and available in your `$PATH`. Place all raw sequence files in your working directory and rename the input files as follows:

- `R1.fastq.gz`
- `R2.fastq.gz`

In the terminal, navigate (`cd`) to your folder containing these files and run:

```
source thesource.txt
```

---

### **Workflow Steps: From FASTQ to Contact Maps**

Below are detailed steps, with the reasoning and logic before each action, followed by the actual command(s) to execute.

---

#### **1. Unique Molecular Identifier (UMI) Extraction**

**Reasoning:**  
Extract the UMI from reads and append it to each sequence's read ID using UMI-tools to allow accurate deduplication later.

**Command:**

```
umi_tools extract -I R1.fastq.gz --bc-pattern=NNNNNNNNNNNNNNNNNNNN --read2-in=R2.fastq.gz --stdout=processed.R1.fastq.gz --read2-out=processed.R2.fastq.gz
```

---

#### **2. Linker Sequence Filtering (Read1/Read2 synchronization via seqkit)**

**Reasoning:**  
Identify reads in R1 that contain the linker sequence, then synchronize and select the corresponding R2 by read ID.

**Commands:**

```
seqkit grep -s -R 23:45 -i -r -p GCTCTTCCGATCT processed.R1.fastq.gz -o linked.processed.R1.fastq.gz
zgrep "@" linked.processed.R1.fastq.gz | sed -n 's/^@\(.*\) .*/\1/p' > linked.rIDs.txt
seqkit grep -f linked.rIDs.txt processed.R2.fastq.gz -o linked.processed.R2.fastq.gz
```

---

#### **3. Seed vs. Receptor Group Assignment (Read2 content-specific separation)**

**Reasoning:**  
Partition reads based on whether Read2 has a specific seed linker sequence, to segregate seed and receptor populations for custom processing.

**Commands:**

```
seqkit grep -s -p CTCATGCGTAGGTAGGCGAC linked.processed.R2.fastq.gz -o linked.processed.R2.seed.fastq.gz
seqkit grep -s -v -p CTCATGCGTAGGTAGGCGAC linked.processed.R2.fastq.gz -o linked.processed.R2.notseed.fastq.gz
```

---

#### **4. Read1 Synchronization Based on Seed/Notseed Assignment**

**Reasoning:**  
Select paired Read1s corresponding to each seed/notseed Read2 by matching read IDs to ensure groups stay in sync.

**Commands:**

```
zgrep "@" linked.processed.R2.seed.fastq.gz | sed -n 's/^@\(.*\) .*/\1/p' > seeds.rIDS.txt
seqkit grep -f seeds.rIDS.txt linked.processed.R1.fastq.gz -o linked.processed.R1.seed.fastq.gz
seqkit grep -v -f seeds.rIDS.txt linked.processed.R1.fastq.gz -o linked.processed.R1.notseed.fastq.gz
```

---

#### **5. Adapter and Quality Trimming**

**Reasoning:**  
Trim reads appropriately based on grouping, optimizing for compatibility with downstream aligners.

**Commands:**

```
trim_galore --fastqc --gzip --clip_R1 40 --paired linked.processed.R1.notseed.fastq.gz linked.processed.R2.notseed.fastq.gz
trim_galore --fastqc --gzip --clip_R1 40 --clip_R2 21 --paired linked.processed.R1.seed.fastq.gz linked.processed.R2.seed.fastq.gz
```

---

#### **6. Alignment and Sorting**

**Reasoning:**  
Align trimmed reads separately for seeds and receptors to the yeast reference genome. Sorting facilitates merging and downstream processing.

**Commands:**

```
bowtie2 -p 16 -x s_cer2011 -1 linked.processed.R1.notseed_val_1.fq.gz -2 linked.processed.R2.notseed_val_2.fq.gz | samtools view -@ 16 -bS - > PE.notseed.bam
samtools sort PE.notseed.bam -o PE.notseed.sorted.bam -@ 16
bowtie2 -p 16 -x s_cer2011 -1 linked.processed.R1.seed_val_1.fq.gz -2 linked.processed.R2.seed_val_2.fq.gz | samtools view -@ 16 -bS - > PE.seed.bam
samtools sort PE.seed.bam -o PE.seed.sorted.bam -@ 16
```

---

#### **7. Merging, Filtering, Resorting, and Indexing**

**Reasoning:**  
Merge aligned reads and filter for high-quality pairs (mapQ ≥30), then re-sort and index for efficient access.

**Commands:**

```
samtools merge merged.bam PE.notseed.sorted.bam PE.seed.sorted.bam -@ 16
samtools view -h -bS -q 30 merged.bam > merged.q.bam -@ 16
samtools sort merged.q.bam -o merged.q.sorted.bam -@ 16
samtools index merged.q.sorted.bam
```

---

#### **8. Deduplication Using UMI-tools**

**Reasoning:**  
Eliminate PCR duplicates by leveraging UMI information.

**Command:**

```
umi_tools dedup -I merged.q.sorted.bam --paired --output-stats --chrom=chr10 -S deduplicated.bam
```

---

#### **9. BAM Indexing for Deduplicated Set**

**Command:**

```
samtools index deduplicated.bam
```

---

#### **10. Read Information Extraction**

**Reasoning:**  
Extract relevant annotation fields for each read using the provided python scripts.

**Command:**

```
python UMI_Reads_bed.py
```

**Output:**  
Generates a `.bed` file:

- Chromosome / start / end / pair length / midpoint / strand / UMI / read ID

---

#### **11. Size Filtering**

**Reasoning:**  
Ensure only biologically relevant read pairs are retained.

**Command:**

```
awk '($4>10 && $4<1000)' read_info.bed > read_info.2.bed
```

---

#### **12. UMI Sorting**

**Command:**

```
sort -t$'\t' -k7,7 read_info.2.bed > read_info.s.bed
```

---

#### **13. Midpoint Delta Calculation and Filtering**

**Reasoning:**  
Filter out potential PCR/sequencing artifacts or ambiguities.

**Commands:**

```
awk ' {print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $5-pv} {pv=$5}' read_info.s.bed > read_info.s.delta.bed
awk '($9>10||$9<0)' read_info.s.delta.bed > read_info.filt.txt
```

---

#### **14. Misalignment Filtering**

**Commands:**

```
awk 'NR==FNR{threshold[$1]=$2; next} $1 in threshold && $5 <= threshold[$1]' sacCer.arabicsize.txt read_info.filt.txt > read_info.filt.thr.txt
awk '$2 >= 0' read_info.filt.thr.txt > read_info.filt.thr2.txt
```

---

#### **15. Pairwise Interaction Generation via Python Script**

**Reasoning:**  
Group by UMI and compute all-vs-all combinations for interaction mapping.

**Command:**

```
python Readinfobed--pairs.py
```

**Output:**  
Whitespace-separated pair file with:

- length1 – chr1 – midpoint1 – mockfrag1 – length2 – chr2 – midpoint2 -mockfrag2

---

#### **16. .hic Map Generation with Juicer Tools**

**Commands:**

```
awk '$2 > $6 {print $5,$6,$7,$8,$1,$2,$3,$4,$9} $2<=$6 {print}' read_info.filt.thr2.txt.juicer | parsort -k2,2d -k6,6d > read_info.s.juicer
java -Xmx16g -jar path-to-your-folder/juicer_tools.jar pre read_info.s.juicer maps.hic sacCer3 -r 10,25,33,50,75,100,150,200,300,400,500,1000,2500,5000,7500,10000,15000,20000,30000,50000,100000,250000,500000,1000000
```

---

#### **17. .mcool Map Generation with Cooler**

**Commands:**

```
awk -F' ' '{print $2 "\t" $3 "\t" $6 "\t" $7 "\t"}'  read_info.filt.thr2.txt.juicer > pairs
bedtools makewindows -g /Users/delamara/juicer/references/sacCer2.chrsize3.txt -w 5 > bins.5.txt
cooler cload pairs -c1 1 -p1 2 -c2 3 -p2 4 bins.5.txt pairs extr.cool
cooler zoomify -r 5N extr.cool
```

---

#### **18. Downstream Analyses (Optional)**

- **Cooltools/coolpup.py**: Generate pile-ups, insulation scores, and distance-decay plots.
- **pyGenomeTracks**: Produce browser tracks.
