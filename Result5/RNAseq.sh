input=/home/liusai/RNA-seq/data/MLT
output=/home/liusai/RNA-seq/result/MLT

mkdir -p $output/{fastp,bam,FC}

cat MLT.txt | while read id
do
    echo "Processing $id"

    ## ---------- fastp ----------
    mkdir -p $output/fastp/${id}
    fastp -w 8 \
      -i $input/${id}.R1.raw.fastq.gz \
      -I $input/${id}.R2.raw.fastq.gz \
      -o $output/fastp/${id}_R1.clean.fq.gz \
      -O $output/fastp/${id}_R2.clean.fq.gz \
      -h $output/fastp/${id}/${id}.html \
      -j $output/fastp/${id}/${id}.json

    ## ---------- HISAT2 + 排序 ----------
    hisat2 --dta -p 32 \
      -x /home/liusai/RNA-sequence/index/hsa \
      -1 $output/fastp/${id}_R1.clean.fq.gz \
      -2 $output/fastp/${id}_R2.clean.fq.gz \
    | samtools sort -@10 -o $output/bam/${id}_hisat2.sorted.bam

    samtools index $output/bam/${id}_hisat2.sorted.bam
done

## ---------- featureCounts ----------
cd $output/FC

featureCounts -T 10 \
-a /home/liusai/RNA-sequence/index/GCF_000001405.40_GRCh38.p14_genomic.gtf \
-o read.count \
-p -B -C \
-t exon -g gene_id \
*_hisat2.sorted.bam
