#!/bin/bash

PROJECT=$1
MDNA_HASH_ADDRESS=$2
SAMPLE_SHEET=$3

REFERENCE="/isilon/sequencing/GATK_resource_bundle/bwa_mem_0.7.5a_ref/human_g1k_v37_decoy.fasta"
P3_1KG="/isilon/sequencing/1000genomes/Full_Project/Sep_2014/20130502/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5.20130502.sites.vcf.gz"
ExAC="/isilon/sequencing/ExAC/Release_0.3/ExAC.r0.3.sites.vep.vcf.gz"

KEY="/isilon/sequencing/CIDRSeqSuiteSoftware/gatk/GATK_2/lee.watkins_jhmi.edu.key"

TIMESTAMP=`date '+%F.%H-%M-%S'`

# set

# REFINE GENOTYPES FOR SNPS USING ALLELE FREQS FROM 1KG AND ExAC

/isilon/sequencing/Kurt/Programs/Java/jdk1.7.0_25/bin/java -jar \
/isilon/sequencing/CIDRSeqSuiteSoftware/gatk/GATK_3/GenomeAnalysisTK-nightly-2015-01-15-g92376d3/GenomeAnalysisTK.jar \
-T CalculateGenotypePosteriors \
-R $REFERENCE \
--variant /isilon/sequencing/Seq_Proj/$PROJECT/JOINT_CALL/$MDNA_HASH_ADDRESS/VQSR/all/variants.all.vcf.gz \
--supporting $P3_1KG \
--supporting $ExAC \
--disable_auto_index_creation_and_locking_when_reading_rods \
-et NO_ET \
-K $KEY \
-o /isilon/sequencing/Seq_Proj/$PROJECT/JOINT_CALL/$MDNA_HASH_ADDRESS/VQSR/variants.all.refined.vcf

# ADD/UPDATE INFO ANNOTATIONS

/isilon/sequencing/Kurt/Programs/Java/jdk1.7.0_25/bin/java -jar \
/isilon/sequencing/CIDRSeqSuiteSoftware/gatk/GATK_3/GenomeAnalysisTK-3.3-0/GenomeAnalysisTK.jar \
-T VariantAnnotator \
-R $REFERENCE \
--variant /isilon/sequencing/Seq_Proj/$PROJECT/JOINT_CALL/$MDNA_HASH_ADDRESS/VQSR/variants.all.refined.vcf \
-L /isilon/sequencing/Seq_Proj/$PROJECT/JOINT_CALL/$MDNA_HASH_ADDRESS/VQSR/variants.all.refined.vcf \
-A SampleList \
-A GCContent \
-A VariantType \
-A GenotypeSummaries \
-A AlleleBalance \
--disable_auto_index_creation_and_locking_when_reading_rods \
-et NO_ET \
-K $KEY \
-o /isilon/sequencing/Seq_Proj/$PROJECT/JOINT_CALL/$MDNA_HASH_ADDRESS/VQSR/$PROJECT".REFINED.vcf"

sed 's/\r//g' $SAMPLE_SHEET \
| awk 'NR>1' \
| cut -d "," -f 1,8,12,18 \
| sort \
| uniq \
| awk 'BEGIN {FS=","} \
{split($2,smtag,"@"); print "qsub","-N","C.01_EXTRACT_VCF_WHOLE_GENOME."smtag[1]"_"smtag[2],\
"-o","/isilon/sequencing/Seq_Proj/"$1"/LOGS/"$2".EXTRACT.VCF.WHOLE.GENOME.MS.log",\
"-e","/isilon/sequencing/Seq_Proj/"$1"/LOGS/"$2".EXTRACT.VCF.WHOLE.GENOME.MS.log",\
"/isilon/sequencing/Seq_Proj/"$1"/SCRIPTS/C.01_EXTRACT_VCF_WHOLE_GENOME.sh",\
$1,$2,$3,$4,"'$MDNA_HASH_ADDRESS'""\n""sleep 3s"}' | bash

sed 's/\r//g' $SAMPLE_SHEET \
| awk 'NR>1' \
| cut -d "," -f 1,8,12,18 \
| sort \
| uniq \
| awk 'BEGIN {FS=","} \
{split($2,smtag,"@"); print "qsub","-N","C.02_EXTRACT_VCF_CODING."smtag[1]"_"smtag[2],\
"-o","/isilon/sequencing/Seq_Proj/"$1"/LOGS/"$2".EXTRACT.VCF.CODING.MS.log",\
"-e","/isilon/sequencing/Seq_Proj/"$1"/LOGS/"$2".EXTRACT.VCF.CODING.MS.log",\
"/isilon/sequencing/Seq_Proj/"$1"/SCRIPTS/C.02_EXTRACT_VCF_CODING.sh",\
$1,$2,$3,$4,"'$MDNA_HASH_ADDRESS'""\n""sleep 3s"}' | bash

sed 's/\r//g' $SAMPLE_SHEET \
| awk 'NR>1' \
| cut -d "," -f 1,8,12,18 \
| sort \
| uniq \
| awk 'BEGIN {FS=","} \
{split($2,smtag,"@"); print "qsub","-N","C.03_EXTRACT_VCF_WHOLE_GENOME_PASS_VARIANT."smtag[1]"_"smtag[2],\
"-o","/isilon/sequencing/Seq_Proj/"$1"/LOGS/"$2".EXTRACT.VCF.WHOLE.GENOME.PASS.VARIANT.MS.log",\
"-e","/isilon/sequencing/Seq_Proj/"$1"/LOGS/"$2".EXTRACT.VCF.WHOLE.GENOME.PASS.VARIANT.MS.log",\
"/isilon/sequencing/Seq_Proj/"$1"/SCRIPTS/C.03_EXTRACT_VCF_WHOLE_GENOME_PASS_VARIANT.sh",\
$1,$2,$3,$4,"'$MDNA_HASH_ADDRESS'""\n""sleep 3s"}' | bash

sed 's/\r//g' $SAMPLE_SHEET \
| awk 'NR>1' \
| cut -d "," -f 1,8,12,18 \
| sort \
| uniq \
| awk 'BEGIN {FS=","} \
{split($2,smtag,"@"); print "qsub","-N","C.04_EXTRACT_VCF_CODING_PASS_VARIANT."smtag[1]"_"smtag[2],\
"-o","/isilon/sequencing/Seq_Proj/"$1"/LOGS/"$2".EXTRACT.VCF.CODING.PASS.VARIANT.MS.log",\
"-e","/isilon/sequencing/Seq_Proj/"$1"/LOGS/"$2".EXTRACT.VCF.CODING.PASS.VARIANT.MS.log",\
"/isilon/sequencing/Seq_Proj/"$1"/SCRIPTS/C.04_EXTRACT_VCF_CODING_PASS_VARIANT.sh",\
$1,$2,$3,$4,"'$MDNA_HASH_ADDRESS'""\n""sleep 3s"}' | bash

sed 's/\r//g' $SAMPLE_SHEET \
| awk 'NR>1' \
| cut -d "," -f 1,8,12,18 \
| sort \
| uniq \
| awk 'BEGIN {FS=","} \
{split($2,smtag,"@"); print "qsub","-N","C.05_EXTRACT_SNV_WHOLE_GENOME."smtag[1]"_"smtag[2],\
"-o","/isilon/sequencing/Seq_Proj/"$1"/LOGS/"$2".EXTRACT.SNV.WHOLE.GENOME.MS.log",\
"-e","/isilon/sequencing/Seq_Proj/"$1"/LOGS/"$2".EXTRACT.SNV.WHOLE.GENOME.MS.log",\
"/isilon/sequencing/Seq_Proj/"$1"/SCRIPTS/C.05_EXTRACT_SNV_WHOLE_GENOME.sh",\
$1,$2,$3,$4,"'$MDNA_HASH_ADDRESS'""\n""sleep 3s"}' | bash

sed 's/\r//g' $SAMPLE_SHEET \
| awk 'NR>1' \
| cut -d "," -f 1,8,12,18 \
| sort \
| uniq \
| awk 'BEGIN {FS=","} \
{split($2,smtag,"@"); print "qsub","-N","C.06_EXTRACT_SNV_CODING."smtag[1]"_"smtag[2],\
"-o","/isilon/sequencing/Seq_Proj/"$1"/LOGS/"$2".EXTRACT.SNV.CODING.MS.log",\
"-e","/isilon/sequencing/Seq_Proj/"$1"/LOGS/"$2".EXTRACT.SNV.CODING.MS.log",\
"/isilon/sequencing/Seq_Proj/"$1"/SCRIPTS/C.06_EXTRACT_SNV_CODING.sh",\
$1,$2,$3,$4,"'$MDNA_HASH_ADDRESS'""\n""sleep 3s"}' | bash

sed 's/\r//g' $SAMPLE_SHEET \
| awk 'NR>1' \
| cut -d "," -f 1,8,12,18 \
| sort \
| uniq \
| awk 'BEGIN {FS=","} \
{split($2,smtag,"@"); print "qsub","-N","C.07_EXTRACT_SNV_WHOLE_GENOME_PASS_VARIANT."smtag[1]"_"smtag[2],\
"-o","/isilon/sequencing/Seq_Proj/"$1"/LOGS/"$2".EXTRACT.SNV.WHOLE.GENOME.PASS.VARIANT.MS.log",\
"-e","/isilon/sequencing/Seq_Proj/"$1"/LOGS/"$2".EXTRACT.SNV.WHOLE.GENOME.PASS.VARIANT.MS.log",\
"/isilon/sequencing/Seq_Proj/"$1"/SCRIPTS/C.07_EXTRACT_SNV_WHOLE_GENOME_PASS_VARIANT.sh",\
$1,$2,$3,$4,"'$MDNA_HASH_ADDRESS'""\n""sleep 3s"}' | bash

sed 's/\r//g' $SAMPLE_SHEET \
| awk 'NR>1' \
| cut -d "," -f 1,8,12,18 \
| sort \
| uniq \
| awk 'BEGIN {FS=","} \
{split($2,smtag,"@"); print "qsub","-N","C.08_EXTRACT_SNV_CODING_PASS_VARIANT."smtag[1]"_"smtag[2],\
"-o","/isilon/sequencing/Seq_Proj/"$1"/LOGS/"$2".EXTRACT.SNV.CODING.PASS.VARIANT.MS.log",\
"-e","/isilon/sequencing/Seq_Proj/"$1"/LOGS/"$2".EXTRACT.SNV.CODING.PASS.VARIANT.MS.log",\
"/isilon/sequencing/Seq_Proj/"$1"/SCRIPTS/C.08_EXTRACT_SNV_CODING_PASS_VARIANT.sh",\
$1,$2,$3,$4,"'$MDNA_HASH_ADDRESS'""\n""sleep 3s"}' | bash

sed 's/\r//g' $SAMPLE_SHEET \
| awk 'NR>1' \
| cut -d "," -f 1,8,12,18 \
| sort \
| uniq \
| awk 'BEGIN {FS=","} \
{split($2,smtag,"@"); print "qsub","-N","C.09_EXTRACT_INDEL_WHOLE_GENOME."smtag[1]"_"smtag[2],\
"-o","/isilon/sequencing/Seq_Proj/"$1"/LOGS/"$2".EXTRACT.INDEL.WHOLE.GENOME.MS.log",\
"-e","/isilon/sequencing/Seq_Proj/"$1"/LOGS/"$2".EXTRACT.INDEL.WHOLE.GENOME.MS.log",\
"/isilon/sequencing/Seq_Proj/"$1"/SCRIPTS/C.09_EXTRACT_INDEL_WHOLE_GENOME.sh",\
$1,$2,$3,$4,"'$MDNA_HASH_ADDRESS'""\n""sleep 3s"}' | bash

sed 's/\r//g' $SAMPLE_SHEET \
| awk 'NR>1' \
| cut -d "," -f 1,8,12,18 \
| sort \
| uniq \
| awk 'BEGIN {FS=","} \
{split($2,smtag,"@"); print "qsub","-N","C.10_EXTRACT_INDEL_CODING."smtag[1]"_"smtag[2],\
"-o","/isilon/sequencing/Seq_Proj/"$1"/LOGS/"$2".EXTRACT.INDEL.CODING.MS.log",\
"-e","/isilon/sequencing/Seq_Proj/"$1"/LOGS/"$2".EXTRACT.INDEL.CODING.MS.log",\
"/isilon/sequencing/Seq_Proj/"$1"/SCRIPTS/C.10_EXTRACT_INDEL_CODING.sh",\
$1,$2,$3,$4,"'$MDNA_HASH_ADDRESS'""\n""sleep 3s"}' | bash

sed 's/\r//g' $SAMPLE_SHEET \
| awk 'NR>1' \
| cut -d "," -f 1,8,12,18 \
| sort \
| uniq \
| awk 'BEGIN {FS=","} \
{split($2,smtag,"@"); print "qsub","-N","C.11_EXTRACT_INDEL_WHOLE_GENOME_PASS_VARIANT."smtag[1]"_"smtag[2],\
"-o","/isilon/sequencing/Seq_Proj/"$1"/LOGS/"$2".EXTRACT.INDEL.WHOLE.GENOME.PASS.VARIANT.MS.log",\
"-e","/isilon/sequencing/Seq_Proj/"$1"/LOGS/"$2".EXTRACT.INDEL.WHOLE.GENOME.PASS.VARIANT.MS.log",\
"/isilon/sequencing/Seq_Proj/"$1"/SCRIPTS/C.11_EXTRACT_INDEL_WHOLE_GENOME_PASS_VARIANT.sh",\
$1,$2,$3,$4,"'$MDNA_HASH_ADDRESS'""\n""sleep 3s"}' | bash

sed 's/\r//g' $SAMPLE_SHEET \
| awk 'NR>1' \
| cut -d "," -f 1,8,12,18 \
| sort \
| uniq \
| awk 'BEGIN {FS=","} \
{split($2,smtag,"@"); print "qsub","-N","C.12_EXTRACT_INDEL_CODING_PASS_VARIANT."smtag[1]"_"smtag[2],\
"-o","/isilon/sequencing/Seq_Proj/"$1"/LOGS/"$2".EXTRACT.INDEL.CODING.PASS.VARIANT.MS.log",\
"-e","/isilon/sequencing/Seq_Proj/"$1"/LOGS/"$2".EXTRACT.INDEL.CODING.PASS.VARIANT.MS.log",\
"/isilon/sequencing/Seq_Proj/"$1"/SCRIPTS/C.12_EXTRACT_INDEL_CODING_PASS_VARIANT.sh",\
$1,$2,$3,$4,"'$MDNA_HASH_ADDRESS'""\n""sleep 3s"}' | bash

sed 's/\r//g' $SAMPLE_SHEET \
| awk 'NR>1' \
| cut -d "," -f 1,8,12,18 \
| sort \
| uniq \
| awk 'BEGIN {FS=","} N==2 \
{split($2,smtag,"@"); print "qsub","-N","C.13_VQSR_SNP_PLOT."smtag[1]"_"smtag[2],\
"-o","/isilon/sequencing/Seq_Proj/"$1"/LOGS/"$2".VQSR.SNP.PLOT.MS.log",\
"-e","/isilon/sequencing/Seq_Proj/"$1"/LOGS/"$2".VQSR.SNP.PLOT.MS.log",\
"/isilon/sequencing/Seq_Proj/"$1"/SCRIPTS/C.13_VQSR_SNP_PLOT_MS.sh",\
$1,$2,$3,$4,"'$MDNA_HASH_ADDRESS'""\n""sleep 3s"}' | bash

sed 's/\r//g' $SAMPLE_SHEET \
| awk 'NR>1' \
| cut -d "," -f 1,8,12,18 \
| sort \
| uniq \
| awk 'BEGIN {FS=","} NR==2 \
{split($2,smtag,"@"); print "qsub","-N","C.14_VQSR_INDEL_PLOT."smtag[1]"_"smtag[2],\
"-o","/isilon/sequencing/Seq_Proj/"$1"/LOGS/"$2".VQSR.INDEL.PLOT.MS.log",\
"-e","/isilon/sequencing/Seq_Proj/"$1"/LOGS/"$2".VQSR.INDEL.PLOT.MS.log",\
"/isilon/sequencing/Seq_Proj/"$1"/SCRIPTS/C.14_VQSR_INDEL_PLOT_MS.sh",\
$1,$2,$3,$4,"'$MDNA_HASH_ADDRESS'""\n""sleep 3s"}' | bash

sed 's/\r//g' $SAMPLE_SHEET \
| awk 'NR>1' \
| cut -d "," -f 1,8,12,18 \
| sort \
| uniq \
| awk 'BEGIN {FS=","} \
{split($2,smtag,"@"); print "qsub","-N","C.15_EXTRACT_VCF_WHOLE_GENOME_PASS_ALL."smtag[1]"_"smtag[2],\
"-o","/isilon/sequencing/Seq_Proj/"$1"/LOGS/"$2".EXTRACT.VCF.WHOLE.GENOME.PASS.ALL.MS.log",\
"-e","/isilon/sequencing/Seq_Proj/"$1"/LOGS/"$2".EXTRACT.VCF.WHOLE.GENOME.PASS.ALL.MS.log",\
"/isilon/sequencing/Seq_Proj/"$1"/SCRIPTS/C.15_EXTRACT_VCF_WHOLE_GENOME_PASS_ALL.sh",\
$1,$2,$3,$4,"'$MDNA_HASH_ADDRESS'""\n""sleep 3s"}' | bash

sed 's/\r//g' $SAMPLE_SHEET \
| awk 'NR>1' \
| cut -d "," -f 1,8,12,18 \
| sort \
| uniq \
| awk 'BEGIN {FS=","} \
{split($2,smtag,"@"); print "qsub","-N","C.16_EXTRACT_VCF_CODING_PASS_ALL."smtag[1]"_"smtag[2],\
"-o","/isilon/sequencing/Seq_Proj/"$1"/LOGS/"$2".EXTRACT.VCF.CODING.PASS.ALL.MS.log",\
"-e","/isilon/sequencing/Seq_Proj/"$1"/LOGS/"$2".EXTRACT.VCF.CODING.PASS.ALL.MS.log",\
"/isilon/sequencing/Seq_Proj/"$1"/SCRIPTS/C.16_EXTRACT_VCF_CODING_PASS_ALL.sh",\
$1,$2,$3,$4,"'$MDNA_HASH_ADDRESS'""\n""sleep 3s"}' | bash

sed 's/\r//g' $SAMPLE_SHEET \
| awk 'NR>1' \
| cut -d "," -f 1,8,12,18 \
| sort \
| uniq \
| awk 'BEGIN {FS=","} \
{split($2,smtag,"@"); print "qsub","-N","C.17_EXTRACT_INDEL_WHOLE_GENOME_PASS_ALL."smtag[1]"_"smtag[2],\
"-o","/isilon/sequencing/Seq_Proj/"$1"/LOGS/"$2".EXTRACT.INDEL.WHOLE.GENOME.PASS.ALL.MS.log",\
"-e","/isilon/sequencing/Seq_Proj/"$1"/LOGS/"$2".EXTRACT.INDEL.WHOLE.GENOME.PASS.ALL.MS.log",\
"/isilon/sequencing/Seq_Proj/"$1"/SCRIPTS/C.17_EXTRACT_INDEL_WHOLE_GENOME_PASS_ALL.sh",\
$1,$2,$3,$4,"'$MDNA_HASH_ADDRESS'""\n""sleep 3s"}' | bash

sed 's/\r//g' $SAMPLE_SHEET \
| awk 'NR>1' \
| cut -d "," -f 1,8,12,18 \
| sort \
| uniq \
| awk 'BEGIN {FS=","} \
{split($2,smtag,"@"); print "qsub","-N","C.18_EXTRACT_INDEL_CODING_PASS_ALL."smtag[1]"_"smtag[2],\
"-o","/isilon/sequencing/Seq_Proj/"$1"/LOGS/"$2".EXTRACT.INDEL.CODING.PASS.ALL.MS.log",\
"-e","/isilon/sequencing/Seq_Proj/"$1"/LOGS/"$2".EXTRACT.INDEL.CODING.PASS.ALL.MS.log",\
"/isilon/sequencing/Seq_Proj/"$1"/SCRIPTS/C.18_EXTRACT_INDEL_CODING_PASS_ALL.sh",\
$1,$2,$3,$4,"'$MDNA_HASH_ADDRESS'""\n""sleep 3s"}' | bash

sed 's/\r//g' $SAMPLE_SHEET \
| awk 'NR>1' \
| cut -d "," -f 1,8,12,18 \
| sort \
| uniq \
| awk 'BEGIN {FS=","} \
{split($2,smtag,"@"); print "qsub","-N","C.19_EXTRACT_SNV_WHOLE_GENOME_PASS_ALL."smtag[1]"_"smtag[2],\
"-o","/isilon/sequencing/Seq_Proj/"$1"/LOGS/"$2".EXTRACT.SNV.WHOLE.GENOME.PASS.ALL.MS.log",\
"-e","/isilon/sequencing/Seq_Proj/"$1"/LOGS/"$2".EXTRACT.SNV.WHOLE.GENOME.PASS.ALL.MS.log",\
"/isilon/sequencing/Seq_Proj/"$1"/SCRIPTS/C.19_EXTRACT_SNV_WHOLE_GENOME_PASS_ALL.sh",\
$1,$2,$3,$4,"'$MDNA_HASH_ADDRESS'""\n""sleep 3s"}' | bash

sed 's/\r//g' $SAMPLE_SHEET \
| awk 'NR>1' \
| cut -d "," -f 1,8,12,18 \
| sort \
| uniq \
| awk 'BEGIN {FS=","} \
{split($2,smtag,"@"); print "qsub","-N","C.20_EXTRACT_SNV_CODING_PASS_ALL."smtag[1]"_"smtag[2],\
"-o","/isilon/sequencing/Seq_Proj/"$1"/LOGS/"$2".EXTRACT.SNV.CODING.PASS.ALL.MS.log",\
"-e","/isilon/sequencing/Seq_Proj/"$1"/LOGS/"$2".EXTRACT.SNV.CODING.PASS.ALL.MS.log",\
"/isilon/sequencing/Seq_Proj/"$1"/SCRIPTS/C.20_EXTRACT_SNV_CODING_PASS_ALL.sh",\
$1,$2,$3,$4,"'$MDNA_HASH_ADDRESS'""\n""sleep 3s"}' | bash

echo
echo START BGZIP-ING NEW VCF FILE AT $TIMESTAMP

bgzip-0.2.6 -c /isilon/sequencing/Seq_Proj/$PROJECT/JOINT_CALL/$MDNA_HASH_ADDRESS/VQSR/$PROJECT".REFINED.vcf" \
>| /isilon/sequencing/Seq_Proj/$PROJECT/JOINT_CALL/$MDNA_HASH_ADDRESS/VQSR/$PROJECT".REFINED.vcf.gz"

tabix-0.2.6 -f -p vcf \
/isilon/sequencing/Seq_Proj/$PROJECT/JOINT_CALL/$MDNA_HASH_ADDRESS/VQSR/$PROJECT".REFINED.vcf.gz"
