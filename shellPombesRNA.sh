#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --gres-flags=enforce-binding
#SBATCH --time=36:00:00
#SBATCH -J sRNAaln
#SBATCH -p cpu-gen8
module load Bowtie
module load FASTX-Toolkit
module load SAMtools
module load BEDTools

#Define directory
raw=/mounts/sarkies/bioc1877/spombesRNA/raw_data/raw_files
ref=/mounts/sarkies/bioc1877/spombesRNA/reference/Schizosaccharomyces_pombe_all_chromosomes.fa


gunzip *.fa.gz
for i in `ls $raw/*.fa`
do
fastx_collapser -i "$i" >col"$i"
perl Fasta_min_length.pl -m 18 col"$i" >tcol_"$i"

bowtie -f -v 0 -k 1 --best -S -x $ref tcol_"$i" > $raw/"$i".sam

samtools view -bS $raw/"$i".sam > $raw/"$i".bam
bamToBed -i $raw/"$i".bam > $raw/"$i".bed
done
