#eddy dowle 2020 pomonella pooled sequencing
#bash code for pooled sequencing:

#clean reads 
#cut adapt/trim
cutadapt \
            -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
            -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
            -o Input_R1.cutadapt.fastq.gz -p Input_R2.cutadapt.fastq.gz \
            Input_R1.fastq.gz Input_R2.fastq.gz
 
java -jar /homes/eddydowle/bin/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 10 -phred33 \
        Input_R1.cutadapt.fastq.gz Input_R1.cutadapt.fastq.gz \
        Input_R1.cutadapt.pairclean.fastq.gz Input_R1.cutadapt.unpairclean.fastq.gz  \
        Input_R2.cutadapt.pairclean.fastq.gz Input_R2.cutadapt.unpairclean.fastq.gz \
        LEADING:10 TRAILING:10 SLIDINGWINDOW:4:15 MINLEN:40


#map to genome
bwa mem 28612_ref_Rhagoletis_zephyria_1.0_chrUn.renamed.fa Input_R1.cutadapt.pairclean.fastq.gz Input_R2.cutadapt.pairclean.fastq.gz -t 6 > Input.cut.trim.bwamem.renamed.sam

#remove secondary/PCR duplicates 

java -jar picard.jar MarkDuplicates \
      I=Input.cut.trim.bwamem.renamed.bam \
      O=Input.cut.trim.bwamem.renamed.bam \
      REMOVE_DUPLICATES=true
      M=marked_dup_metrics.txt

#estimate coverage
bedtools genomecov -ibam Input.cut.trim.bwamem.renamed.bam -g 28612_ref_Rhagoletis_zephyria_1.0_chrUn.renamed.fa > Input.coverage.txt
 
samtools depth Input.cut.trim.bwamem.renamed.bam | awk '{sum+=$3;cnt++}END{print sum/cnt" "sum}' > Input_samtoolsdepth
 
java -jar /homes/eddydowle/bin/picard-tools-2.2.1/picard.jar CollectAlignmentSummaryMetrics \
      R=28612_ref_Rhagoletis_zephyria_1.0_chrUn.renamed.fa \
      I=Input.cut.trim.bwamem.renamed.bam \
      O=Input.picardTools

#realign around your indels
#first find your indels in a merged bam file of all 6 samples
#used bamtools for this as GATK painfully slow - bamtools shouldnt be used to call indels as it cant deal with poolseq but it will find the indels ok
#at this point we only need to be able to find the indels so we can realign the reads around them

#run in array
# 'grep' to grab the chromosome names out of the fasta file remove the '>' and then use 'split' to break it up to 3000 lines per file
#then call indels/snps (cant just do indels) on each chromosome individually via a array job
#here Im using a SGE server for slurm so its $SLURM_ARRAY_TASK_ID and set #SBATCH --array=1-29 #to number of total jobs in the header

mkdir mpileup.$SGE_TASK_ID
for line in $(<28612_ref_Rhagoletis_zephyria_1.0_chrUn.renamed.chr.$SGE_TASK_ID); do samtools mpileup -d 250 -m 1 -E --min-MQ 30 --min-BQ 30 --BCF --output-tags AD,DP,DV,SP -r "${line}" -f 28612_ref_Rhagoletis_zephyria_1.0_chrUn.renamed.fa -o mpileup.$SGE_TASK_ID/"${line}_Merged.6files.cut.trim.bwamem.renamed.sort.nosecondary.remdup.rg.picard.sort.bcf" Merged.6files.cut.trim.bwamem.renamed.sort.nosecondary.remdup.rg.picard.sort.bam
done

#this creates a lot of files in a lot of folders
#merge them up
#merge within each folder
OUTPUT="$(tr '\n' ' ' < 28612_ref_Rhagoletis_zephyria_1.0_chrUn.renamed.chr.bcf.$SGE_TASK_ID)"
bcftools concat $OUTPUT -o 28612_ref_Rhagoletis_zephyria_1.0_chrUn.renamed.chr.bcf.$SGE_TASK_ID.merged.bcf

#merge across folders (no longer in a array job single job flick the flag of!)
OUTPUT="$(tr '\n' ' ' < final.bcffiles)"
bcftools concat $OUTPUT -o 28612_ref_Rhagoletis_zephyria_1.0_chrUn.renamed.chr.bcf.AllChromosome.merged.vcf

#index and create a file with JUST variable indel calls:
bcftools index 28612_ref_Rhagoletis_zephyria_1.0_chrUn.renamed.chr.bcf.AllChromosome.merged.vcf
bcftools call --skip-variants snps --multiallelic-caller --variants-only  -O v -o 28612_ref_Rhagoletis_zephyria_1.0_chrUn.renamed.chr.bcf.AllChromosome.merged.indel.vcf 28612_ref_Rhagoletis_zephyria_1.0_chrUn.renamed.chr.bcf.AllChromosome.merged.vcf

#realign your bam files with gatk using the indel calls:
java -jar picard-tools-2.2.1/picard.jar SortVcf \
       I=28612_ref_Rhagoletis_zephyria_1.0_chrUn.renamed.chr.bcf.AllChromosome.merged.vcf \
       O=28612_ref_Rhagoletis_zephyria_1.0_chrUn.renamed.chr.bcf.AllChromosome.merged.sort.vcf

#I arrayed this across the bams again
for line in $(<bamfileList.$SGE_TASK_ID); do  java -jar GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar \
   -T IndelRealigner \
   -R 28612_ref_Rhagoletis_zephyria_1.0_chrUn.renamed.fa \
   -I "${line}.bam" \
   -known 28612_ref_Rhagoletis_zephyria_1.0_chrUn.renamed.chr.bcf.AllChromosome.merged.indel.sort.vcf \
   -targetIntervals 28612_ref_Rhagoletis_zephyria_1.0_chrUn.renamed.chr.bcf.AllChromosome.merged.indel.sort.intervals \
   -o "${line}.indelrealign.bam"
done

#now its secondary/pcr duplicate free and realinged around the indels
#ANGSD is generally run for genotype probabilities, but we want just the allele counts at each base and a MAF based of allele count
#ANGSD works by having a text file with the path and file name of the a bam file for each of your samples (they cant be in the same bam)
#everything ANGSD will output will be in the order of your bam file list.

#first find sites in angsd
#so first lets just find the sites

angsd -b bamfilesPoolSeqlist -nThreads 6 -GL 2 \
        -doMajorMinor 2 -minQ 30 -minMapQ 30 -minInd 6 -doVcf 1 -minMaf 0.05 \
        -setMinDepthInd 10 -nind 6 -doCounts 1 \
        -dumpCounts 4 -domaf 8 -doGeno 32 -doPost 1 -doGlf \
        -out PoolSeqAllSamplesDepthInd10Qual30minind6.maf005.domajmin1

#-GL 2 (do genotype likelihoods using GATK approach)
#-doMajorMinor 2 (use counts to find Maj/Min)
#-minQ (miniumum site quality)
#-minMapQ (miniumum map quality) 
#-minInd 6 (just set it to your number of individuals ~ you could set it to less than but then you will have to deal with missing data, Im not sure if its worth the hassel on the poolseq)
#-doVcf 1 (will output a vcf file)
#-minMaf 0.05 
#-setMinDepthInd 10 (ten counts per indidividual at a site)
#-nind 6 (total individuals)
#-doCounts 1 (count reads)
#-dumpCounts 4 (output the counts)
#-doMaf 8 (calculate the allele frequency of the counts not genotype likelihood)
#-doGeno 32 (output genotype likelihoods for PCA)
#-doPost 1 
#-doGlf (will output a beagle file)
#-out (outfile ~there will be multiple files tagged to this)

#this will give you 
# XX.arg (just a log file)
# XX.counts.gz (counts file)
# XX.mafs.gz (allele frequency estiamtes file)
# XX.geno.gz (genotype likelihoods file)
# XX.pos.gz
# XX.vcf.gz

#ok so you want to 

zcat XX.mafs.gz > XX.mafs | awk '{print $1"\t"$2"\t"$3"\t"$4}' > sites.file

#this is the sites angsd found open it up and remove the first line and have a look

#it should look something like this:
head sites.file
NW_016157085.1	1000119	G	C
NW_016157085.1	100013	A	T
NW_016157085.1	1000139	G	C
NW_016157085.1	1000144	G	C
NW_016157085.1	1000157	A	G
NW_016157085.1	1000175	T	C
NW_016157085.1	1000199	A	T
NW_016157085.1	1000228	T	C
NW_016157085.1	1000230	A	G
NW_016157085.1	1000252	G	C

#we are going to feed this back into angsd
#this makes things much faster etc, angsd no longer has to check each site it just goes through the list
#it also means we can set the major minor as the same in all the six pools

#first up index the sites file and create a file of the 'chromosomes' (scaffolds really)
cut -f1 sites.file |sort|uniq >sites.file.chrs
angsd sites index sites.file

#now lets re-run angsd and we are now going to set things up for poolseq

#for each pool run (example here AppleAve):
angsd -b bamfileAppleAve -nThreads 6 -GL 2 \
        -doMajorMinor 3 -minQ 30 -minMapQ 30 -minInd 1 \
        -setMinDepthInd 10 -nind 1 -doCounts 1 \
        -dumpCounts 3 -domaf 8  -doPost 1 -doVcf 1 \
        -out PoolSeqAppleAveDepthInd10Qual30.sites -sites sites.file -rf sites.file.chrs

#remember -doMajorMinor 3 incase the maj/min are around the other way

#then the .mafs.gz file has the allele frequency estimate per site for each of the pools

#ok now we have all our files we need to sort things out
#we want to do a fisher test on the count data but we need a file that looks like:

#chr	pos	majlist	minlist	appleavemaj	appleavemin	hawavemaj	hawavemin
#NW_016157085.1	117	A	T	17	0	23	0
#NW_016157085.1	127	A	G	8	10	17	8
#NW_016157085.1	131	T	A	9	8	13	13
#NW_016157085.1	133	A	G	10	8	15	13
#NW_016157085.1	139	G	T	19	0	27	0
#NW_016157085.1	145	T	G	14	7	22	5
#NW_016157085.1	154	T	A	18	4	28	2
#NW_016157085.1	164	T	C	14	9	23	7
#NW_016157085.1	174	A	C	18	4	24	5

#I wrote a script called poolseq.rawtomajmin.py 
#That will take the poolseq output and build something we can use:

#chromo	position	major	minor	AppleAve_A	AppleAve_C	AppleAve_G	AppleAve_T	AppleEarly_A	AppleEarly_C	AppleEarly_G	AppleEarly_T	AppleLate_A	AppleLate_C	AppleLate_G	AppleLate_T	HawAve_A	HawAve_C	HawAve_G	HawAve_T	HawEarly_A	HawEarly_C	HawEarly_G	HawEarly_T	HawLate_A	HawLate_C	HawLate_G	HawLate_T
#NW_016157085.1	117	A	T	17	0	0	0	19	0	0	2	15	0	0	3	23	0	0	0	22	0	0	0	28	0	0	4
#NW_016157085.1	127	A	G	8	0	10	0	16	0	3	0	12	0	4	0	17	0	8	0	19	0	5	0	30	0	2	0
#NW_016157085.1	131	T	A	8	0	0	9	6	0	0	16	5	1	0	11	13	0	0	13	12	0	0	15	16	0	0	17
#NW_016157085.1	133	A	G	10	0	8	0	14	0	8	0	10	0	8	0	15	0	13	0	14	0	11	0	13	0	20	0
#NW_016157085.1	139	G	T	0	0	19	0	1	0	24	2	2	0	12	4	2	0	27	0	0	0	31	0	0	0	33	3
#NW_016157085.1	145	T	G	0	0	7	14	0	0	3	24	0	0	4	16	0	0	5	22	0	0	7	24	0	0	3	32
#NW_016157085.1	154	T	A	4	0	0	18	3	0	0	29	4	0	0	21	2	0	0	28	2	0	0	34	2	0	0	30
#NW_016157085.1	164	T	C	0	9	0	14	0	9	0	26	0	4	0	23	0	7	0	23	0	10	0	28	0	8	0	28
#NW_016157085.1	174	A	C	18	4	0	0	36	5	0	0	19	2	0	0	24	5	0	0	25	12	0	0	27	7	0	0

#chr	pos	majlist	minlist	appleavemaj	appleavemin	appleaveTriple	appleearlymaj	appleearlymin	appleearlyTriple	applelatemaj	applelatemin	applelateTriple	hawavemaj	hawavemin	hawaveTriple	hawearlymaj	hawearlymin	hawearlyTriple	hawlatemaj	hawlatemin	hawlateTriple
#NW_016157085.1	117	A	T	17	0	None	19	2	None	15	3	None	23	0	None	22	0	None	28	4	None
#NW_016157085.1	127	A	G	8	10	None	16	3	None	12	4	None	17	8	None	19	5	None	30	2	None
#NW_016157085.1	131	T	A	9	8	None	16	6	None	11	5	None	13	13	None	15	12	None	17	16	None
#NW_016157085.1	133	A	G	10	8	None	14	8	None	10	8	None	15	13	None	14	11	None	13	20	None
#NW_016157085.1	139	G	T	19	0	None	24	2	None	12	4	None	27	0	None	31	0	None	33	3	None
#NW_016157085.1	145	T	G	14	7	None	24	3	None	16	4	None	22	5	None	24	7	None	32	3	None
#NW_016157085.1	154	T	A	18	4	None	29	3	None	21	4	None	28	2	None	34	2	None	30	2	None
#NW_016157085.1	164	T	C	14	9	None	26	9	None	23	4	None	23	7	None	28	10	None	28	8	None
#NW_016157085.1	174	A	C	18	4	None	36	5	None	19	2	None	24	5	None	25	12	None	27	7	None

#using awk build a file like this:
#chr	pos	majlist	minlist	appleavemaj	appleavemin	hawavemaj	hawavemin
#NW_016157085.1	117	A	T	17	0	23	0
#NW_016157085.1	127	A	G	8	10	17	8
#NW_016157085.1	131	T	A	9	8	13	13
#NW_016157085.1	133	A	G	10	8	15	13
#NW_016157085.1	139	G	T	19	0	27	0
#NW_016157085.1	145	T	G	14	7	22	5
#NW_016157085.1	154	T	A	18	4	28	2
#NW_016157085.1	164	T	C	14	9	23	7
#NW_016157085.1	174	A	C	18	4	24	5

#then just run the other script in that folder on it FisherExactTestPoolmajminslot.py 
#this one just runs a fisher test on the counts at each line

#annotation:
#the rhago genome is up on a software called snpEff 
#this is why you need the vcf we arnt going to use the likelihood scores but snpEff needs the site information in a vcfformat
#snpeff will annotate whether our sites are missense, intronic etc 
java -Xmx200G -jar snpEff.jar Rhagoletis_zephyria PoolSeqAllSamplesDepthInd10Qual30minind6.domajmin3.sites.vcf > PoolSeqAllSamplesDepthInd10Qual30minind6.domajmin3.sites.vcf

#you will need to use vcfSift to make this readable

#data was then built into mysql database for analyses


