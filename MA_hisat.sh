#!/bin/bash -l
#SBATCH --job-name=hisat_featurecounts
#SBATCH --output=output_%j.txt
#SBATCH --error=errors_%j.txt
#SBATCH --time=08:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=2G
#SBATCH --array=0-8
#SBATCH -p cpu-gen8
#SBATCH --mail-type=END
#SBATCH --mail-user=mariana.villalbadelapena@bioch.ox.ac.uk

module load all/Miniconda3/25.7.0-2
conda init
conda activate featurecounts
module load all/SAMtools/1.18-GCC-12.3.0


#Define reference genome
referenceP=Schizosaccharomyces_pombe_all_chromosomes.fa
#Define the output folder
output=hisat_output

input=hisat_output
output2=featurecounts_output

#all_samples 
samples=(Anc L1_T10 L1_T5 L2_T10 L2_T5 L6_T10 L6_T5 L7_T10 L7_T5)

sample=${samples[$SLURM_ARRAY_TASK_ID]} #Note that bash arrays are 0-index based

### Index de genomes (only has to be done once)
#IDXp=$referenceP
#IDXc=$referenceC
 
#hisat2-build $IDXp $IDXp
#hisat2-build $IDXc $IDXc

#samtools faidx $IDXp
#samtools faidx $IDXc
#########################################

hisat2 -x $referenceP -1 ${sample}_1.fq.gz -2 ${sample}_2.fq.gz | samtools sort > $output/${sample}.bam

#sense feature counts
featureCounts -p --countReadPairs -a Schizosaccharomyces_pombe_all_chromosomes.gtf -s 1 -o $output2/${sample}_sensepombecounts.txt $input/${sample}.bam 
featureCounts -p --countReadPairs -a Schizosaccharomyces_pombe_all_chromosomes.gtf  -o $output2/${sample}_pombecounts.txt $input/${sample}.bam
