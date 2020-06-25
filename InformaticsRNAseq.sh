#eddy dowle
#pomonella RNA sequencing

#clean remove adapters
for line in $(<fastqfilesRhago); 
do 
java -jar trimmomatic-0.36.jar SE -phred33 "${line}_R1_001.fastq.gz" "${line}_R1_001.clean.fastq.gz"  ILLUMINACLIP:TruSeq3.fa:2:30:10 LEADING:10 TRAILING:10 SLIDINGWINDOW:4:15 MINLEN:50
done

#star mapping/rsem:
STAR --runMode genomeGenerate --genomeDir genomeSTAR --genomeFastaFiles 28612_ref_Rhagoletis_zephyria_1.0_chrUn.fa --sjdbGTFfile ref_Rhagoletis_zephyria_1.0_top_level.gtf
for line in $(<fastqfilesRhago); _R1_001.clean.fastq.gz
do 
zcat "${line}R1_001.clean.fastq.gz" > "${line}R1_001.clean.fastq"
STAR --genomeDir genomeSTAR --readFilesIn "${line}R1_001.clean.fastq" --runThreadN 14 --outFileNamePrefix "${line}.pubgenome.bam" --sjdbGTFtagExonParentTranscript ref_Rhagoletis_zephyria_1.0_top_level.gff3 --quantMode TranscriptomeSAM --outSAMtype BAM SortedByCoordinate 
rsem-calculate-expression -p 12 --bam --no-bam-output "${line}.pubgenome.bam" 28612_ref_Rhagoletis_zephyria_1.0_chrUn.renamed "${line}.rsem"
done

