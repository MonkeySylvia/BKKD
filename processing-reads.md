Processing the reads
================

``` r
#Add novogene adaptor (didnt figure out the index thing)
GATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG
```

**Trim adaptors**

``` r
for sample in `cat sample.txt`; do
    echo "${sample}_1.fq" "${sample}_2.fq"
    java -jar /programs/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 30 \
    -phred33 "${sample}_1.fq" "${sample}_2.fq" \
    "${sample}_1_clean_paired_qc.fastq" "${sample}_1_clean_unpaired_qc.fastq" \
    "${sample}_2_clean_paired_qc.fastq" "${sample}_2_clean_unpaired_qc.fastq" \
    ILLUMINACLIP:/workdir/nc499/BKKD/process/Trimmomatic-0.39_combined_adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
done
```

**map**

``` r
source /programs/HISAT2/hisat2.sh
#alignment
for sample in `cat sample.txt`; do
  hisat2 -q --max-intronlen 750000 --rna-strandness FR --fr -p 30 -k 100 \
  -x /workdir/nc499/reference/danRer11.nonalt.ht2 \
  -1 ${sample}_1_clean_paired_qc.fastq \
  -2 ${sample}_2_clean_paired_qc.fastq \
  -S ${sample}.sam;
  done > hisat.log &
```

**TE count** <br> *need python 2, so I added python2 in TEcount and changed it to TEcount2*

``` r
for sample in `ls *.sam`; do
  TEcount2 -b $sample \
  --GTF /workdir/nc499/reference/dr11_ref/Danio_rerio.GRCz11.98.chr.gtf \
  --TE /workdir/nc499/BKKD/danRer11_rmsk_TE.gtf \
  --format SAM --mode multi \
  --project $sample \
  --stranded reverse; done
```

``` r
sed -i 's/\.sam//' *.sam.cntTable

paste Exp1_8.sam.cntTable <(cut -f2 Exp1_12.sam.cntTable) <(cut -f2 Exp2_12.sam.cntTable) <(cut -f2 Exp3_8.sam.cntTable) <(cut -f2 Sc1_8.sam.cntTable) <(cut -f2 Sc1_12.sam.cntTable) <(cut -f2 Sc2_12.sam.cntTable) <(cut -f2 Sc3_8.sam.cntTable) | sed 's/ /\t/g' > 8samples.cntTable
```
