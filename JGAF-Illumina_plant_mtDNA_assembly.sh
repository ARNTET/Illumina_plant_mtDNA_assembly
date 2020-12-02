#!/bin/bash
#SBATCH --partition=fast
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=2000
#SBATCH --job-name="Illumina_Assembly"
#SBATCH --output=JGAF-Illumina_Assembly%j.out

# Modules needed ###############################################################
module load trimmomatic
module load fastq-join
module load velvet
module load SPAdes
module load samtools

# Sequences needed #############################################################
ref_mtDNA=LACWENDEL_mtDNA
bwa index $ref_mtDNA".fasta"
ILLUMINACLIP=NexteraPE-PE.fa

# Sample #######################################################################

#A# TRIMMING
mkdir a_trimmed
trimmomatic \
PE \
$sample"_R1.fastq.gz" \
$sample"_R2.fastq.gz" \
a_trimmed/$sample".R1.fq.gz" \
a_trimmed/$sample".R1U.fq.gz" \
a_trimmed/$sample".R2.fq.gz" \
a_trimmed/$sample".R2U.fq.gz" \
ILLUMINACLIP:NexteraPE-PE.fa:2:30:10 SLIDINGWINDOW:4:15 MINLEN:75

##B## MERGING
mkdir b_joined
fastq-join \
a_trimmed/$sample".R1.fq.gz" \
a_trimmed/$sample".R2.fq.gz" \
-o b_joined/$sample".R1.fq" \
-o b_joined/$sample".R2.fq" \
-o b_joined/$sample".joined.fq"

#C# 1st ASSEMBLY
velveth \
c_assemblyV \
93 \
-shortPaired \
-fastq \
-separate \
b_joined/$sample".R1.fq" \
b_joined/$sample".R2.fq"
velvetg \
c_assemblyV \
-scaffolding no

#D# MAPPING CONTIGS TO REFERENCE GENOME
mkdir d_mapped_contigs
bwa mem -t 4 \
$ref_mtDNA".fasta"\
c_assemblyV/contigs.fa \
> d_mapped_contigs/mapped_contigs.sam
samtools view -b \
d_mapped_contigs/mapped_contigs.sam \
> d_mapped_contigs/mapped_contigs.bam
samtools fasta -F4 \
d_mapped_contigs/mapped_contigs.bam \
> d_mapped_contigs/mapped_contigs.fasta


##F## ASSEMBLY
mkdir -p e_assemblyS
spades.py \
-k 33,55,77,99,127 \
-t 16 \
--cov-cutoff auto \
--trusted-contigs d_mapped_contigs/mapped_contigs.fasta \
-1 b_joined/$sample".R1.fastq" \
-2 b_joined/$sample".R2.fastq" \
-s b_joined/$sample".joined.fastq" \
-o e_assemblyS/$sample".spades"

##F## FILTERING
# => Contigs need to be blasted to identify mitochondrial coverage
mkdir -p f_filtering
samtools faidx e_assemblyS/$sample".spades"/contigs.fasta
wc -l e_assembly/$sample".spades"/contigs.fasta.fai

# X = Upper limit to mitochondrial coverage
# Y = Lower limit to mitochondrial coverage

cat e_assembly/$sample".spades"/contigs.fasta.fai \
| sed -e $'s/_/\t/g' \
| sed 's/\./,/g' \
| awk '$6 < X && $6 > Y' \
| cut -f 1-6 \
| sed -e $'s/\t/_/g' \
| sed 's/,/\./g' > f_filtering/$sample"list_contigs.txt"

xargs samtools faidx \
e_assemblyS/$sample".spades"/contigs.fast \
< f_filtering/$sample"list_contigs.txt" \
> f_filtering/$sample"contigs_flt.fasta"

##G## SCAFFOLDING
mkdir -p g_scaffolding
cd g_scaffolding

perl ~/sspace_basic/SSPACE_Basic.pl \
-l libraries.txt \
-s f_filtering/$sample"contigs_flt.fasta" \
-x 0 -m 32 -o 20 -t 0 -k 5 -a 0.70 -n 15 -p 0 -v 0 -z 0 -g 0 -T 1 \
-b $sample".scaffolds.fasta"
