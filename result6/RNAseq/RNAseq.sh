input=/home/liusai/RNA-seq/data/MLT
output=/home/liusai/RNA-seq/result/MLT
mkdir -p $output/fastp
mkdir -p $output/sam
mkdir -p $output/bam
mkdir -p $output/FC

cat MLT.txt | while read id
do
    mkdir -p $output/fastp/${id}
    fastp -i $input/${id}.R1.raw.fastq.gz -o $output/fastp/${id}_filter_R1.fq.gz  \
          -I $input/${id}.R2.raw.fastq.gz -O $output/fastp/${id}_filter_R2.fq.gz \
          -h $output/fastp/${id}/${id}_report.html

    mkdir -p $output/sam/${id}
    hisat2 -t -p 32 -x /home/liusai/RNA-sequence/index/hsa \
           -1 $output/fastp/${id}_filter_R1.fq.gz \
           -2 $output/fastp/${id}_filter_R2.fq.gz \
           -S $output/sam/${id}/${id}.sam

    samtools view -@10 -bS $output/sam/${id}/${id}.sam > $output/FC/${id}.bam

done
featureCounts -T 10 -a /home/liusai/RNA-sequence/index/GCF_000001405.40_GRCh38.p14_genomic.gtf -o read.count -p -t exon -g gene_id  *.bam


