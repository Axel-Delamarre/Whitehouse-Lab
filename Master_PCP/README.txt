

### PCP data processing workflow

Requirements
- UMI-tools https://umi-tools.readthedocs.io/en/latest/
- seqkit https://bioinf.shenwei.me/seqkit/
- trim_galore https://github.com/FelixKrueger/TrimGalore
- bowtie2
- samtools
- python
- juicer_tools https://github.com/aidenlab/JuicerTools
- cooler https://github.com/open2c/cooler
- bedtools

Visualize the .hic fiels in Juicerbox https://github.com/aidenlab/Juicebox
Visualize the .mcool file in HiGlass https://higlass.io/ or using pyGenomeTracks https://pygenometracks.readthedocs.io/en/latest/



This prpocess can be run through source:

- Copy all these files in you folder, rename the .fastq.gz into R1.fastq.gz and R2.fastq.gz
- In the terminal cd to the appropriate folder and run
source thesource.txt



Steps in the data processing from fastq files to .hic and .mcool maps.

•	The UMI is extracted and appended to the read ID using UMI-tools 
umi_tools extract -I R1.fastq.gz --bc-pattern=NNNNNNNNNNNNNNNNNNNN --read2-in=R2. fastq.gz --stdout=processed.R1.fastq.gz --read2-out=processed.R2. fastq.gz
•	The reads containing the linker sequences on Read1 are selected using seqkit, the corresponding Read2 are selected
    seqkit grep -s -R 23:45 -i -r -p GCTCTTCCGATCT processed.R1.fastq.gz -o linked.processed.R1.fastq.gz
    zgrep "@" linked.processed.R1.fastq.gz | sed -n 's/^@\(.*\) .*/\1/p' > linked.rIDs.txt
    seqkit grep -f linked.rIDs.txt processed.R2.fastq.gz -o linked.processed.R2.fastq.gz
•	The reads are separated in seed or receptor group based on the presence of the seed linker in their sequence on Read2:
    seqkit grep -s -p CTCATGCGTAGGTAGGCGAC linked.processed.R2.fastq.gz -o linked.processed.R2.seed.fastq.gz
    seqkit grep -s -v -p CTCATGCGTAGGTAGGCGAC linked.processed.R2.fastq.gz -o linked.processed.R2.notseed.fastq.gz
•	The corresponding R1 are processed through read ID identification:
    zgrep "@" linked.processed.R2.seed.fastq.gz | sed -n 's/^@\(.*\) .*/\1/p' > seeds.rIDS.txt
    seqkit grep -f seeds.rIDS.txt linked.processed.R1.fastq.gz -o linked.processed.R1.seed.fastq.gz
    seqkit grep -v -f seeds.rIDS.txt linked.processed.R1.fastq.gz -o linked.processed.R1.notseed.fastq.gz
•	The reads are trimmed differently on Read2 if they are seed or not:
    trim_galore --fastqc --gzip --clip_R1 40 --paired linked.processed.R1.notseed.fastq.gz linked.processed.R2.notseed.fastq.gz
    trim_galore --fastqc --gzip --clip_R1 40 --clip_R2 21 --paired linked.processed.R1.seed.fastq.gz linked.processed.R2.seed.fastq.gz
•	The seeds and receptors are aligned separately and sorted on the yeast genome 
    bowtie2 -p 16 -x s_cer2011 -1 linked.processed.R1.notseed_val_1.fq.gz -2 linked.processed.R2.notseed_val_2.fq.gz | samtools view -@ 16 -bS - > PE.notseed.bam
    samtools sort PE.notseed.bam -o PE.notseed.sorted.bam -@ 16
    bowtie2 -p 16 -x s_cer2011 -1 linked.processed.R1.seed_val_1.fq.gz -2 linked.processed.R2.seed_val_2.fq.gz | samtools view -@ 16 -bS - > PE.seed.bam
    samtools sort PE.seed.bam -o PE.seed.sorted.bam -@ 16
•	The reads are merged, quality filtered, resorted and indexed
    samtools merge merged.bam PE.notseed.sorted.bam PE.seed.sorted.bam -@ 16
    samtools view -h -bS -q 30 merged.bam > merged.q.bam -@ 16
    samtools sort merged.q.bam -o merged.q.sorted.bam -@ 16
    samtools index merged.q.sorted.bam
•	The reads are deduplicated using UMI-tools
    umi_tools dedup -I merged.q.sorted.bam --paired --output-stats --chrom=chr10 -S deduplicated.bam
•	The bam file is re indexed
    samtools index deduplicated.bam
•	The read information are extracted through a python script that will be available on github
    python UMI_Reads_bed.py
      - This produces a bed file with the following columns:
	    Chromosome / start / end / pair length / midpoint / strand / UMI / read ID
•	Reads are selecting for their length to be comprised between 10bp and 1Kb
    awk '($4>10 && $4<1000)' read_info.bed > read_info.2.bed
•	Reads are sorted on the UMI sequence
    sort -t$'\t' -k7,7 read_info.2.bed > read_info.s.bed
•	The delta in midpoint of following reads is calculated, the reads for which the midpoint is too closed are discarded. This is important for 2 reasons (1- UMI-tools might have miss duplicates. 2- In the PCP product, each nucleosome can produce 2 molecules emanating from its 2 strands, only one of them must be used in the analysis).
    awk ' {print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $5-pv} {pv=$5}' read_info.s.bed > read_info.s.delta.bed
    awk '($9>10||$9<0)' read_info.s.delta.bed > read_info.filt.txt
•	Reads are filtered for any misalignements
    awk 'NR==FNR{threshold[$1]=$2; next} $1 in threshold && $5 <= threshold[$1]' sacCer.arabicsize.txt read_info.filt.txt > read_info.filt.thr.txt
    awk '$2 >= 0' read_info.filt.thr.txt > read_info.filt.thr2.txt
•	Reads are grouped by UMI sequenced and paiwise interactions combination are deduced through a python script that will be available on github. The midpoint of the read is used as a value. Importantly, if the molecule is bigger than 260bp, the value used is the dyad of the nucleosomes closest to the tag extremity: Tag extremity +/- (depending on the strand) 83 (half a nucleosome)
    python Readinfobed--pairs.py
      - This produce a pair file, white-space separated,  with the following columns
	      length1 – chr1 – midpoint1 – mockfrag1 – length2 – chr2 – midpoint2 -mockfrag2
•	The produced file is used to generate .hic maps using juicertools pre:
    awk '$2 > $6 {print $5,$6,$7,$8,$1,$2,$3,$4,$9} $2<=$6 {print}' read_info.filt.thr2.txt.juicer | parsort -k2,2d -k6,6d > read_info.s.juicer
    java -Xmx16g -jar path-to-your-folder/juicer_tools.jar pre read_info.s.juicer maps.hic sacCer3 -r 10,25,33,50,75,100,150,200,300,400,500,1000,2500,5000,7500,10000,15000,20000,30000,50000,100000,250000,500000,1000000
•	The produced files is used to generate .mcool maps using cooler:
    awk -F' ' '{print $2 "\t" $3 "\t" $6 "\t" $7 "\t"}'  read_info.filt.thr2.txt.juicer > pairs 
    bedtools makewindows -g /Users/delamara/juicer/references/sacCer2.chrsize3.txt -w 5 > bins.5.txt
    cooler cload pairs -c1 1 -p1 2 -c2 3 -p2 4 bins.5.txt pairs extr.cool 
    cooler zoomify -r 5N extr.cool

•	Cooltools and coolpup.py were used to generates Pile-up analysis, Insulation score, Distance-decay plots.
•	pyGenomeTracks was used to generates the genome browser tracks. 
