# mtb_ExpEvo

Scripts for data analyzing sequence data from experimental evolution associated with the following publication: eLife2022;11:e78454 (DOI: https://doi.org/10.7554/eLife.78454)


# popoolationStats.py
This script takes in alignment files (BAMs) from pool-seq, subsamples to minimize effects of differential coverage and calculates genome wide population statistics (pi and Tajima's D) using Popoolation.

# popoolationSyncToTSV.py
Converts Popoolation2 sync file to TSV with allele frequencies. Filters mutations based on bed file, allele count, coverage and frequency across samples.

# slidingWindowCoverage.py
This script calculates coverage in sliding windows from BAM(s) and outputs compiled information.
