patha=/home/liusai/Cut-tag/data/MLT/clean_data
pathr=/home/liusai/Cut-tag/RESULT/MLT
path1=$patha
path2=$pathr/1filterdata
path3=$pathr/2bam
path4=$pathr/3bed
path5=$pathr/4macs2ou
path7=$pathr/7bw

#mkdir $pathr
mkdir $path2
mkdir $path3
mkdir $path4
mkdir $path5
mkdir $path7
cat $patha/MLT.txt | while read id
  do
      mkdir $path2/${id}
      mkdir $path3/${id}
      mkdir $path5/${id}
      fastp -i $path1/${id}.clean.1.fastq.gz -o $path2/${id}/${id}_filter_R1.fq.gz  -I  $path1/${id}.clean.2.fastq.gz -O $path2/${id}/${id}_filter_R2.fq.gz -h $path2/${id}/${id}_report.html
      fastqc -o $path2/${id} -t 6 $path2/${id}/${id}_filter_R1.fq.gz $path2/${id}/${id}_filter_R2.fq.gz
      bowtie2 --end-to-end --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 -p 10\
              -x /home/liusai/index/bowtie2/human/ncbi/p14\
              -1 $path2/${id}/${id}_filter_R1.fq.gz\
              -2 $path2/${id}/${id}_filter_R2.fq.gz\
              -S ${path3}/${id}/${id}_bowtie2.sam
   samtools view -@10 -bS -F 0x04 ${path3}/${id}/${id}_bowtie2.sam > ${path3}/${id}/${id}_bowtie2.bam
   sambamba markdup -t 6 -r -p ${path3}/${id}/${id}_bowtie2.bam ${path3}/${id}/${id}_bowtie2f.bam
   samtools sort -@ 8 ${path3}/${id}/${id}_bowtie2f.bam -o ${path3}/${id}/${id}_bowtie2sort.bam
   samtools index -@ 8 ${path3}/${id}/${id}_bowtie2sort.bam
   macs2 callpeak -t  ${path3}/${id}/${id}_bowtie2sort.bam -p 1e-5 -f BAMPE -g hs -n ${id} -B --outdir ${path5} > ${path3}/${id}/${id}.txt
   conda activate deeptools
   bamCoverage -b  ${path3}/${id}/${id}_bowtie2sort.bam -of bigwig -o ${path7}/${id}.bw -p 20 --ignoreDuplicates --binSize 10 --normalizeUsing RPKM
   conda activate rna
done
