# Assembly of plant mitochondrial DNA using a Illumina MiSeq paired-end reads and a reference genome

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

In this study Illumina paired-end short reads sequencing strategy was used to build the mtDNA of 9 accessions from Lactuca spp.. To do so, total leaf DNA of WT and radA plants was quantified with a QuBit Fluorometer (Life Technologies) and libraries were prepared with the Nextera XT Library kit, according to manufacturer’s recommendations (Illumina) using 100 ng of each DNA sample. Final libraries were quantified, checked on a Bioanalyzer 2100 (Agilent) and sequenced on an Illumina Miseq system (2 × 150 paired-end reads).

### Mitochondrial assembly

We built our pipeline based on a recent work from Garcia et al., 2019 (see JGAF-Illumina_plant_mtDNA_assembly.sh). Illumina reads from Lactuca species were trimmed with Trimmomatic to remove adapter sequences and low quality ends and joined with Fastq-join (github.com/ExpressionAnalysis/ea-utils) when paired reads mapped to each other. Then those reads were assembled into contigs with a first assembler: Velvet. Those contigs were mapped on the L. sativa var. capitata L. nidus tenerrima mtDNA sequence and mapped Velvet contigs were used as trusted contigs for a second round of assembly with SPAdes. SPAdes contigs were identified manually by BLAST and filtered according to their coverage. Filtered SPAdes contigs were extended with SSPACE into long scaffolds (up 170 kb). Finally, contigs were manually assembled into single circular molecules. 
