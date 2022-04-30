#genome indexing
STAR --runThreadN 8 --runMode genomeGenerate --genomeDir /mnt/z/EGA_DATA/Genome_data/GRCh38_STARindex --genomeFastaFiles /mnt/z/EGA_DATA/Genome_data/GRCh38.d1.vd1.fa --sjdbGTFfile /mnt/z/EGA_DATA/Genome_data/gencode.v22.annotation.gtf
#sample handling
PATH_origin="/mnt/z/EGA_DATA/EGAD00001000223/"
cd $PATH_origin
for i in $(ls *_1.rnaseq.fastq.gz)
do
cd $PATH_origin
SAMPLE=$(basename $i "_1.rnaseq.fastq.gz")
EXT=".rnaseq.fastq.gz"
FOR="_1"
REV="_2"
STAR_INDEX="/mnt/z/EGA_DATA/Genome_data/GRCh38_STARindex"
PATH_OUT_COUNT="/mnt/z/EGA_DATA/Output_count/"
PATH_WORK="/mnt/z/EGA_DATA/Work/"
cp $SAMPLE$FOR$EXT $PATH_WORK"$SAMPLE$FOR$EXT"
cp $SAMPLE$REV$EXT $PATH_WORK"$SAMPLE$REV$EXT"
cd $PATH_WORK
echo $SAMPLE.trimming
#trimmomatic
trimmomatic PE -threads 8  $SAMPLE$FOR$EXT $SAMPLE$REV$EXT $SAMPLE$FOR.P.fq $SAMPLE$FOR.U.fq $SAMPLE$REV.P.fq $SAMPLE$REV.U.fq LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20
###STAR
echo $SAMPLE.STAR
STAR  --readFilesIn $SAMPLE$FOR.P.fq $SAMPLE$REV.P.fq --alignIntronMax 1000000 --alignIntronMin 20 --alignMatesGapMax 1000000 --alignSJDBoverhangMin 1 --alignSJoverhangMin 8  --alignSoftClipAtReferenceEnds Yes --chimJunctionOverhangMin 15  --chimMainSegmentMultNmax 1  --chimOutType Junctions SeparateSAMold WithinBAM SoftClip  --chimSegmentMin 15  --genomeDir $STAR_INDEX   --genomeLoad NoSharedMemory --limitSjdbInsertNsj 1200000 --outFileNamePrefix $SAMPLE --outFilterIntronMotifs None  --outFilterMatchNminOverLread 0.33  --outFilterMismatchNmax 999  --outFilterMismatchNoverLmax 0.1  --outFilterMultimapNmax 20  --outFilterScoreMinOverLread 0.33 --outFilterType BySJout  --outSAMattributes NH HI AS nM NM ch  --outSAMstrandField intronMotif  --outSAMtype BAM Unsorted  --outSAMunmapped Within --quantMode TranscriptomeSAM GeneCounts --runThreadN 8 --twopassMode Basic --sjdbGTFfile /mnt/z/EGA_DATA/Genome_data/gencode.v22.annotation.gtf
####HTSEQ
echo $SAMPLE.htseq-count
htseq-count -f bam -r name -s no -a 10 -t exon -i gene_id -m intersection-nonempty $SAMPLE"Aligned.out.bam" /mnt/z/EGA_DATA/Genome_data/gencode.v22.annotation.gtf > $SAMPLE.count.txt
cp $SAMPLE.count.txt $PATH_OUT_COUNT$SAMPLE.count.txt
rm -r *
done
