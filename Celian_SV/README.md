# 302 candidate regions: PCA, Relate, and DoLoReS evidence

This directory contains the compact set of files needed to inspect the 302 candidate regions and the shared-selection hypothesis between European farmed and North American farmed salmon.

## Files

- `merged_302_regions_window_PCA.pdf`: Local-PCA plots for the 302 regions and their analysed windows.
- `302_loci_dotplots.pdf`: dotplots for the 302 candidate regions across the long-read assemblies and outgroups.
- `merged_302_regions_window_haplotype_heatmaps.pdf`: haplotype heatmaps for the 302 regions and their analysed windows.
- `302_region_evidence_summary.tsv`: one-row-per-region evidence summary. It contains the region ID and `chrom/start/end`, Local-PCA detection fields for the combined, farmed, and wild European/North American groups, Relate signal counts and positions for European farmed, European wild, North American farmed, and North American wild, and DoLoReS signal counts and intervals.
- `LocalPCA_302_regions_analysis_qc.tsv`: region-level SNP and PCA-window QC used to account for region size and SNP observation opportunity.
- `Relate_DetectSelection_significant_p1e-4.tsv`: significant Relate positions with population labels.
- `DoLoReS_significant_p1e-4.tsv`: significant DoLoReS intervals and signal strengths.
- `Relate_shared_SNPs_European_farmed_vs_North_American_farmed_302.tsv`: exact shared Relate SNP positions between European farmed and North American farmed for each region, including the number of signal positions in each population, the number of identical positions, their coordinates, and the shared fraction.

The central comparison is whether European farmed and North American farmed show Relate signal at the same positions within a region. Raw signal counts should be interpreted together with region width, filtered SNP count, and analysed window count.
