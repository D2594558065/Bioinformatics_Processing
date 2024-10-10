###构建index
###建立索引
STAR --runMode genomeGenerate \
--runThreadN 10 \
--genomeDir /home/mnt1/wanghao/project/STAR_NMR_index \
--genomeFastaFiles /home/mnt1/wanghao/project/Naked_Mole-Rat_genmoe/GCF_000247695.1_HetGla_female_1.0_genomic.fna \
--sjdbGTFfile /home/mnt1/wanghao/project/Naked_Mole-Rat_genmoe/GCF_000247695.1_HetGla_female_1.0_genomic.gtf \
--sjdbOverhang 99 \
--limitGenomeGenerateRAM  406188524128

####去接头

for i in `cat sample.txt`;
do trim_galore -j 10 -q 15   --stringency 3 --length 20 -e 0.1  --paired ${i}_R1_001.fastq.gz  ${i}_R2_001.fastq.gz  --gzip -o ../trim_glore --fastqc; done

###比对

do STAR \
    --genomeDir /home/mnt1/wanghao/project/STAR_NMR_index \
     --runThreadN 10 \
    --readFilesIn ${i}_R1_001_val_1.fq.gz ${i}_R2_001_val_2.fq.gz  \
    --readFilesCommand zcat \
    --outFileNamePrefix /home/mnt1/wanghao/project/help/CYcgas/fruit_fly/test/NMR/$i \
    --outSAMtype BAM SortedByCoordinate \
    --outBAMsortingThreadN 10 ;
done

###计数
featureCounts -p --countReadPairs -T 10  -t exon -g gene_id -a /home/mnt1/wanghao/project/Naked_Mole-Rat_genmoe/GCF_000247695.1_HetGla_female_1.0_genomic.gtf  -o ./matrix/combined.count ./*.out.bam
