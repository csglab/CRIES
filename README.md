# CRIES
## Counting Reads for Intronic and Exonic Segments

This pipeline describes the process used by the CSG lab for counting RNA-seq reads that map to intronic and exonic segments, primarily for the purpose of measuring mRNA stability.

## Requirements
- Unix-compatible OS
- [R](http://www.r-project.org/) version 3.0.1 or later
- [HISAT2](https://ccb.jhu.edu/software/hisat2/index.shtml)
- [HTSeq](http://www-huber.embl.de/users/anders/HTSeq/doc/overview.html)
- [DESeq](http://bioconductor.org/packages/release/bioc/html/DESeq.html) library
- [REMBRANDTS](https://github.com/csglab/REMBRANDTS)

## Step 1. Creating GTF annotation files

We usually extract the coordinates of intronic and exonic segments from the latest Ensembl annotation file. For the human genome, release 87 of Ensembl annotations can be found [here](ftp://ftp.ensembl.org/pub/release-87/gtf/homo_sapiens) (Homo_sapiens.GRCh38.87.chr.gtf.gz). For the purpose of inferring stability from RNA-seq, we only consider constitutive exons, i.e. exons that are present in all isoforms of a gene. We also consider only isoforms that are supported by both the Ensembl and Havana annotations, as provided in the GTF file.

Unzip the annotation file, and then use a shell script such as this to obtain the coordinates of intronic and constitutive exonic segments of the genes:

```bash
# This is the ource annotation file
# Note that the annotation file should be in Ensembl format (downloaded from Ensembl FTP)
# Specifically, the chromosome names should lack the "chr" prefix
# Also, the first tag in column 9 must be gene_id, and the third tag must be transcript_id
# Also, the source annotation file must have a transcript_source tage in column 9

input="Homo_sapiens.GRCh38.87.chr.gtf"

# 1. Identify transcripts that are supported by both Ensembl and Havana annotations
# 2. Identify and write constitutive exons, i.e. those that appear in all Ensembl/Havana isoforms of a gene
output="./Homo_sapiens.GRCh38.87.chr.consExons.gtf"
cat $input | grep 'transcript_source "ensembl_havana"' | awk -v FS='\t' '$3=="exon" { exonName=$1":"$4":"$5":"$7; split($9, fields, ";"); geneName=fields[1]; transcriptName=fields[3]; printf("%s\t%s\t%s\n",exonName,geneName,transcriptName); }' | sort | uniq | awk -v FS='\t' '{ eCount[$1]++; tCount[$3]++; exonHost[$1]=$2; if(tCount[$3]==1) gCount[$2]++; } END { for(i in eCount) if(eCount[i]==gCount[exonHost[i]]) { split(i,fields,":"); printf("%s\tensembl_havana\texon\t%s\t%s\t.\t%s\t.\t%s;\n",fields[1],fields[2],fields[3],fields[4],exonHost[i]); } }' > $output

# 1. Identify transcripts that are supported by both Ensembl and Havana annotations
# 2. Find all exons
# 3. Write the coordinate of intronic regions, i.e. the regions that separate two adjacent exons of the same gene
output="./Homo_sapiens.GRCh38.87.chr.Introns.gtf"
cat $input | grep 'transcript_source "ensembl_havana"' | awk -v FS='\t' '$3=="exon" { exonName=$1":"$4":"$5":"$7; split($9, fields, ";"); geneName=fields[1]; transcriptName=fields[3]; printf("%s\t%s\t%s\n",exonName,geneName,transcriptName); }' | sort | uniq | awk -v FS='\t' '{ eCount[$1]++; tCount[$3]++; exonHost[$1]=$2; if(tCount[$3]==1) gCount[$2]++; } END { for(i in eCount) { split(i,fields,":"); printf("%s\tensembl_havana\texon\t%s\t%s\t.\t%s\t.\t%s;\n",fields[1],fields[2],fields[3],fields[4],exonHost[i]); } }' | bedtools sort -i stdin | awk -v FS='\t' '{ if( last_exon[$9]==1 && (last_exon_end[$9]+1)<($4-1) ) printf("%s\t%s\tintron\t%i\t%i\t%s\t%s\t%s\t%s\n",$1,$2,last_exon_end[$9]+1,$4-1,$6,$7,$8,$9); last_exon[$9]=1; last_exon_end[$9]=$5; }' > $output
```

## Step 2. Mapping reads

We use HISAT2 for mapping RNA-seq reads, keeping only alignments with MAPQ score ≥30.

```bash
hisat2 -t -p 16 -x <indexPath> -U <fastqFiles> | samtools sort -T <scratchPath> -o <sampleName>.srt.bam -
samtools view -b -F 1548 -q 30 -o <sampleName>.srt.filtered.bam <sampleName>.srt.bam
rm <sampleName>.srt.bam
```

Note that the sort command is designed for single-end sequencing data. For paired-end reads, use option `-n`.

## Step 3. Counting reads that map to intronic or exonic segments of each gene

We use HTSeq-count for counting reads. For counting exonic reads, we run the HTSeq-count using the "intersection-strict" mode, to ensure that the reads that are counted entirely fit within the mature mRNA sequence, and do not overlap alternatively spliced exons. For intronic reads, we run the HTSeq-count using the "union" mode, since any read that even partially overlaps an intronic region likely originates from pre-mature RNA.

```bash
htseq-count -m intersection-strict -f bam -s <strand> <sampleName>.srt.filtered.bam ./Homo_sapiens.GRCh38.87.chr.consExons.gtf > <sampleName>".consExons.counts.txt"
htseq-count -m union -f bam -s <strand> <sampleName>.srt.filtered.bam ./Homo_sapiens.GRCh38.87.chr.Introns.gtf > <sampleName>".Introns.counts.txt"
```

The parameter `<strand>` depends on the strandedness of the library preparation kit. For Illumina TruSeq RNA Library Prep Kit, the correct parameter is often `reverse`. If not sure, try all the three different options `no`, `yes` and `reverse`, and decide accordingly.
Note that for paired-end reads, `-r name` must be used along with BAM files that are sorted by read name.

## Step 4. Normalization

We generally use our package [REMBRANDTS](https://github.com/csglab/REMBRANDTS) for read count normalization, gene filtering, and removing bias in order to obtain stability measurements. REMBRANDTS internally uses DESeq, and performs linear regression of Δexon–Δintron vs. Δintron to obtain residuals that correspond to unbiased estimates of stability.

In case you do not want to use REMBRANDTS, you can use variance-stabilized transformation of DESeq or DESeq2 for read count normalization (for each set of intronic and exonic reads separately), and then take Δexon–Δintron between any two samples as the measure of differential stability. However, note that this measure can be biased, overestimating the stability of transcriptionally down-regulated genes and underestimating the stability of transcriptionally up-regulated genes.
