#!/bin/bash
#SBATCH --partition=fast
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=2000

# Modules needed ###############################################################
module load trimmomatic
module load fastq-join
module load velvet
module load SPAdes
module load samtools

# /!\ To fill /!\ ##############################################################
samples=()
lengthS=${#samples[@]}
ref_mtDNA=
adapters=

# Indexing reference ###########################################################
bwa index $ref_mtDNA

# Step one #######################################################################
a=0
while (($a<$lengthS)); do
  #A# Trimming
  mkdir a_trimmed
  trimmomatic \
  PE \
  ${samples[a]}(a)"_R1.fastq.gz" \
  ${samples[a]}"_R2.fastq.gz" \
  a_trimmed/${samples[a]}".R1.fq.gz" \
  a_trimmed/${samples[a]}".R1U.fq.gz" \
  a_trimmed/${samples[a]}".R2.fq.gz" \
  a_trimmed/${samples[a]}".R2U.fq.gz" \
  ILLUMINACLIP:$adapters:2:30:10 SLIDINGWINDOW:4:15 MINLEN:75
  #B# Merging overlaping reads
  mkdir b_joined
  fastq-join \
  a_trimmed/${samples[a]}".R1.fq.gz" \
  a_trimmed/${samples[a]}".R2.fq.gz" \
  -o b_joined/${samples[a]}".R1.fq" \
  -o b_joined/${samples[a]}".R2.fq" \
  -o b_joined/${samples[a]}".joined.fq"
  #C# First assembly
  velveth \
  c_assemblyV \
  93 \
  -shortPaired \
  -fastq \
  -separate \
  b_joined/${samples[a]}".R1.fq" \
  b_joined/${samples[a]}".R2.fq"
  velvetg \
  c_assemblyV \
  -scaffolding no
  #D# Mapping and filtering mitochondrial contigs
  mkdir d_mapped_contigs
  bwa mem -t 16 \
  $ref_mtDNA \
  c_assemblyV/contigs.fa \
  > d_mapped_contigs/${samples[a]}"mapped_contigs.sam"
  samtools view -b \
  d_mapped_contigs/${samples[a]}"mapped_contigs.sam" \
  > d_mapped_contigs/${samples[a]}"mapped_contigs.bam"
  samtools fasta -F4 \
  d_mapped_contigs/${samples[a]}"mapped_contigs.bam" \
  > d_mapped_contigs/${samples[a]}"mapped_contigs.fasta"
  #E# Second assembly
  mkdir e_assemblyS
  cd e_assemblyS
  spades.py \
  -k 33,55,77,99,127 \
  -t 16 \
  --cov-cutoff auto \
  --trusted-contigs ../d_mapped_contigs/${samples[a]}"mapped_contigs.fasta" \
  -1 ../b_joined/${samples[a]}".R1.fastq" \
  -2 ../b_joined/${samples[a]}".R2.fastq" \
  -s ../b_joined/${samples[a]}".joined.fastq" \
  -o ${samples[a]}".spades"
  cd ..
  samtools faidx e_assemblyS/${samples[a]}".spades"/contigs.fasta
  wc -l e_assembly/${samples[a]}".spades"/contigs.fasta.fai
  let "a=a+1"
done
echo Now contigs need to be blasted to identify mitochondrial coverage

# /!\ To fill /!\ ##############################################################
a=0
while (($a<$lengthS)); do
  echo What is the upper limit to mitochondrial coverage in ${samples[a]}
  read ${samples[a]}"_Upper"
  echo What is the lower limit to mitochondrial coverage in ${samples[a]}
  read ${samples[a]}"_Lower"
  let "a=a+1"
done
echo Thank you for these informations

# Step two #######################################################################
a=0
while (($a<$lengthS)); do
  #F# Filtering
  mkdir f_filtering
  cat e_assembly/${samples[a]}".spades"/contigs.fasta.fai \
  | sed -e $'s/_/\t/g' \
  | sed 's/\./,/g' \
  | awk '$6 < '${samples[a]}"_Upper"' && $6 > '${samples[a]}"_Lower"'' \
  | cut -f 1-6 \
  | sed -e $'s/\t/_/g' \
  | sed 's/,/\./g' > f_filtering/${samples[a]}"list_contigs.txt"
  xargs samtools faidx \
  e_assemblyS/${samples[a]}".spades"/contigs.fast \
  < f_filtering/${samples[a]}"list_contigs.txt" \
  > f_filtering/${samples[a]}"contigs_flt.fasta"
  #G# Scaffolding
  mkdir g_scaffolding
  cd g_scaffolding
  perl SSPACE_Basic.pl \
  -l libraries.txt \
  -s f_filtering/${samples[a]}"contigs_flt.fasta" \
  -x 0 -m 32 -o 20 -t 0 -k 5 -a 0.70 -n 15 -p 0 -v 0 -z 0 -g 0 -T 1 \
  -b ${samples[a]}".sspace"
  cd ..
  let "a=a+1"
done
echo DING Your mitochondrial scaffolds are now available 
