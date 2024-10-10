#######去接头
for i in `cat sample.txt`;
do trim_galore -j 4 -q 20  --stringency 3 --length 20 -e 0.1  --paired ${i}_R1_001.fastq.gz  ${i}_R2_001.fastq.gz  --gzip -o ../trimgloreData;done

#######单独处理7号样本 因为文件名称不一样
trim_galore -j 4 -q 20  --stringency 3 --length 20 -e 0.1  --paired WMT-240727-7_S0_L000_R1_000.fastq.gz  WMT-240727-7_S0_L000_R2_000.fastq.gz  --gzip -o ../trimgloreData
#######比对
STAR \
    --genomeDir  /home/mnt1/wanghao/project/index/STAR/CRcm38_mouse \
     --runThreadN 10 \
    --readFilesIn WMT-240727-7_S0_L000_R1_000_val_1.fq.gz WMT-240727-7_S0_L000_R2_000_val_2.fq.gz  \
    --readFilesCommand zcat \
    --outFileNamePrefix ../star_output/WMT-240727-7_S0_L000 \
    --outSAMtype BAM SortedByCoordinate \
    --outBAMsortingThreadN 10 ;


################构建索引
STAR --runMode genomeGenerate \
--runThreadN 10 \
--genomeDir /home/mnt1/wanghao/project/index/STAR/CRcm38_mouse \
--genomeFastaFiles /home/mnt1/wanghao/project/GRCm38/GRCm38.primary_assembly.genome.fa \
--sjdbGTFfile /home/mnt1/wanghao/project/gtf/gencode.vM10.annotation.gtf \
--sjdbOverhang 99 \
--limitGenomeGenerateRAM  406188524128


####进行比对
for i in `cat sample.txt`;
do STAR \
    --genomeDir  /home/mnt1/wanghao/project/index/STAR/CRcm38_mouse \
     --runThreadN 10 \
    --readFilesIn ${i}_R1_001_val_1.fq.gz ${i}_R2_001_val_2.fq.gz  \
    --readFilesCommand zcat \
    --outFileNamePrefix ../star_output/$i \
    --outSAMtype BAM SortedByCoordinate \
    --outBAMsortingThreadN 10 ;
done

#####进行计数
mkdir /home/mnt1/wanghao/project/LINE1_evolution/OR_receptor/LINE1_OR_Data/bulk-RNA/skin/trimglore/matrix
###计数
featureCounts -p -T 10  -t exon -g gene_id -a /home/mnt1/wanghao/project/gtf/gencode.vM10.annotation.gtf -o ../matrix/combined.count  *.out.bam
