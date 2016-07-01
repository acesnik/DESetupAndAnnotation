for BAM in *.sorted.bam
do
        echo $BAM >> sorted_bam_flagstat.txt
        samtools flagstat $BAM >> sorted_bam_flagstat.txt
done