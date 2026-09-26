#!/bin/bash -l
#SBATCH --job-name=featurecounts
#SBATCH --output=output_%j.txt
#SBATCH --error=errors_%j.txt
#SBATCH --time=08:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=2G
#SBATCH -p cpu-gen8
#SBATCH --mail-type=END
#SBATCH --mail-user=mariana.villalbadelapena@bioch.ox.ac.uk


input=hisat_output
output2=featurecounts_output

#featureCounts -p -a Schizosaccharomyces_pombe_all_chromosomes.gtf -o $output2/pombe_counts.txt \
#    $input/L41_pombe.bam $input/L42_pombe.bam $input/L43_pombe.bam $input/L44_pombe.bam \
#    $input/L45_pombe.bam $input/L46_pombe.bam $input/L47_pombe.bam $input/L48_pombe.bam \
#    $input/L49_pombe.bam $input/L50_pombe.bam $input/L51_pombe.bam $input/L52_pombe.bam \
#    $input/L53_pombe.bam $input/L54_pombe.bam $input/L55_pombe.bam $input/L56_pombe.bam \
#    $input/L57_pombe.bam $input/L58_pombe.bam $input/L59_pombe.bam $input/L60_pombe.bam

#featureCounts -p -a Caenorhabditis_elegans.WBcel235.115.gtf -o $output2/elegans_counts.txt \
#     $input/L41_elegans.bam $input/L42_elegans.bam $input/L43_elegans.bam $input/L44_elegans.bam \
#     $input/L45_elegans.bam $input/L46_elegans.bam $input/L47_elegans.bam $input/L48_elegans.bam \
#     $input/L49_elegans.bam $input/L50_elegans.bam $input/L51_elegans.bam $input/L52_elegans.bam \
#     $input/L53_elegans.bam $input/L54_elegans.bam $input/L55_elegans.bam $input/L56_elegans.bam \
#     $input/L57_elegans.bam $input/L58_elegans.bam $input/L59_elegans.bam $input/L60_elegans.bam




featureCounts -p -a Schizosaccharomyces_pombe_all_chromosomes.gtf -o $output2/pombe_countsYES.txt \
    $input/L61_pombe.bam $input/L62_pombe.bam $input/L63_pombe.bam $input/L64_pombe.bam \
    $input/L65_pombe.bam $input/L66_pombe.bam $input/L67_pombe.bam $input/L68_pombe.bam \
    $input/L69_pombe.bam $input/L70_pombe.bam $input/L71_pombe.bam $input/L72_pombe.bam \
    $input/L73_pombe.bam $input/L74_pombe.bam $input/L75_pombe.bam $input/L76_pombe.bam \
    $input/L77_pombe.bam $input/L78_pombe.bam $input/L79_pombe.bam $input/L80_pombe.bam

featureCounts -p -a Caenorhabditis_elegans.WBcel235.115.gtf -o $output2/elegans_countsYES.txt \
     $input/L61_elegans.bam $input/L62_elegans.bam $input/L63_elegans.bam $input/L64_elegans.bam \
     $input/L65_elegans.bam $input/L66_elegans.bam $input/L67_elegans.bam $input/L68_elegans.bam \
     $input/L69_elegans.bam $input/L70_elegans.bam $input/L71_elegans.bam $input/L72_elegans.bam \
     $input/L73_elegans.bam $input/L74_elegans.bam $input/L75_elegans.bam $input/L76_elegans.bam \
     $input/L77_elegans.bam $input/L78_elegans.bam $input/L79_elegans.bam $input/L80_elegans.bam


