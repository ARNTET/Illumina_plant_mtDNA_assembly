# Illumina_plant_mtDNA_assembly
Mitochondrial genome assembly using Illumina MiSeq paired-end reads

# Assembly

## Materials

### Reference genome used in this study:
Lactuca sativa var. capitata L. nidus tenerrima: 

### Raw reads fastq files are available on ebi.ac.uk/ena under the accession numbers:

- Lactuca sativa:
  - UC12100: ERS5267071
  - LAC004500: ERS5267070

- Lactuca saligna:
  - LAC008020:	ERS5267072

- Lactuca serriola:
  - LAC005780: ERS5267075
  - CGN004799: ERS5267073

- Lactuca virosa:
  - LAC006941: ERS5267078
  - CGN013357: ERS5267076
  - CGN019045: ERS5267077



## Methods

### Illumina sequencing

In this study Illumina paired-end short reads sequencing strategy was used to determine the impact of the Arabidopsis thaliana radA mutant on the cpDNA. To do so, total leaf DNA of WT and radA plants was quantified with a QuBit Fluorometer (Life Technologies) and libraries were prepared with the Nextera Flex Library kit, according to manufacturer’s recommendations (Illumina) using 100 ng of each DNA sample. Final libraries were quantified, checked on a Bioanalyzer 2100 (Agilent) and sequenced on an Illumina Miseq system (2 × 150 paired-end reads).

### Coverage analysis

For Illumina sequence analysis of the cpDNA, reads were aligned against the Arabidopsis reference genomes using BWA and filtered to keep only those mapping to the cpDNA. For analysis of rearranged sequences, reads properly pairing to the cpDNA were filtered to only keep those showing a short-clipping sequence (threshold 20 nucleotides) and no indel (looking for the presence of an S in the CIGAR string without any I, D and H). The short-clipping sequences were then extracted with SE-MEI/extractSoftclipped (github.com/dpryan79/SE-MEI) and aligned using bowtie2 against the cpDNA (see JGAF-Atha_cpDNA_SoftClipped.sh). The positions of the short-clipping sequences mapping the cpDNA and of their relative read were rounded down to the upper kb to analyze the localization of the rearrangement (see JGAF-Atha_cpDNA_SoftClipped.R). Those corresponding to the isomerization that results from the recombination involving the large inverted repeats were filtered out.
