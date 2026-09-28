## Repository files
- Epimutations
- **`FASTA_min_length.pl`**: Perl script that filters FASTA sequences by minimum length and updates sequence headers with the first base and sequence length.
- **`shellPombeRNA.sh`**: Bash/SLURM script for processing and aligning small RNA sequencing reads to the *S. pombe* genome using FASTX-Toolkit, Bowtie, SAMtools, and BEDTools.
- **`filter_bed.R`**: R script for retaining 22–26 nt siRNAs and normalizing siRNA data using DESeq.
- **`CUT_OFF.R`**: R script for calculating the Z-score threshold based on technical replicates.
- **`non-parametric_epimutations.R`**: R script for calling epimutations, calculating divergence, and determining the number of epimutations and transitions.
- **`survival.R`**: R script for estimating the mid-life and stability of epimutations using the `survival` and `survminer` packages.
- **`genomic_distribution.R`**: R script for analyzing the genomic distribution of siRNAs and epimutations using the `GenomicDistributions` package.
- **`cluster_profiler.R`**: R script for performing over-representation analysis (ORA) on epimutations using `clusterProfiler`.
- RNA_seq
- **`clusterProfiler_mRNA.R`**: R script for performing over-representation analysis (ORA) on mRNA data using `clusterProfiler`.
- **`MA_hisat.sh`**: Shell script for RNA-seq alignment using HISAT.
- **`feature_counts.sh`**: Shell script for generating read counts from HISAT-aligned RNA-seq data.
- **`mutants_RNAseq.R`**: R script for RNA-seq analysis of mutant strains.
- Annotation files
- **`org.Spombe.egBASE.Rd`**, **`org.Spombe.egORGANISM.Rd`**, **`orgSpombe.eg_dbconn.Rd`**: Files required to create the *S. pombe* annotation database used for enrichment analysis.
- Growth_rate
- **`mutants_growthcurve.R`**: R script for analyzing differences in growth parameters between mutant and wild-type strains using data from `gc_mutantout_edit.csv`.

Raw sequencing files are available in NCBI under the under the accession number: PRJNA1535449
