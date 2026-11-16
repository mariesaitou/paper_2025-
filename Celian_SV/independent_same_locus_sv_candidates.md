# Candidates for independent SV-like changes at the same locus

This note records 302-region candidates where the dotplot suggests that salmon assemblies carry different structural configurations within the same Local-PCA region. The search uses the existing 302-region summaries and the available per-region dotplot and haplotype-heatmap images. The criterion is a candidate for follow-up: the figures can show assembly-specific breaks, contig joins, or orientations, while independent origin requires breakpoint-level comparison.

## PCA_SV_ID_0203 — ssa16:33,928,303–73,490,590

This is the strongest current candidate. The region contains 327 overlapping BayesTyper SVs and 32 direction-discordant comparisons, with a concordance fraction of 0.4483. The dotplot shows a complex multi-segment architecture across the salmon assemblies, including assembly-specific breaks and changes around the GARN-1/GARN-2 configurations. Relate contains 1,530 signals and DoLoReS contains 4 signals. The pattern is consistent with an old structural background plus later lineage-specific rearrangements or recurrent changes at the same broad locus.

- Dotplot: [PCA_SV_ID_0203](scenario_figures/dotplot/PCA_SV_ID_0203.png)
- Haplotype heatmap: [PCA_SV_ID_0203](scenario_figures/heatmap/PCA_SV_ID_0203.png)
- Local-PCA and integrated metrics: [302_region_integrated_metrics.tsv](302_region_integrated_metrics.tsv)

## PCA_SV_ID_0283 — ssa27:10,112,401–16,192,809

This is a high-priority candidate because the salmon assemblies separate into visibly different alignment architectures. European assemblies contain mixed or separated segments, while the North American assemblies show a different predominantly same-orientation configuration. The region contains 49 overlapping BayesTyper SVs and 7 direction-discordant comparisons, with a concordance fraction of 0.3636. Local-PCA cluster separation is 16.8768; Relate contains 443 signals and DoLoReS contains 1 signal. This supports lineage-associated structural configurations at one locus, subject to breakpoint matching.

- Dotplot: [PCA_SV_ID_0283](scenario_figures/dotplot/PCA_SV_ID_0283.png)
- Haplotype heatmap: [PCA_SV_ID_0283](scenario_figures/heatmap/PCA_SV_ID_0283.png)
- Local-PCA and integrated metrics: [302_region_integrated_metrics.tsv](302_region_integrated_metrics.tsv)

## PCA_SV_ID_0276 — ssa26:14,630,742–27,656,297

This region shows a shared broad alignment background with assembly-specific shifted segments, breaks, and alternative orientations. It contains 125 overlapping BayesTyper SVs and 9 direction-discordant comparisons, with a concordance fraction of 0.5000. Relate contains 968 signals and DoLoReS contains 1 signal. The appropriate interpretation is a shared ancestral structural background with possible lineage-specific modifications; breakpoint-level comparison is required to test recurrent independent origin.

- Dotplot: [PCA_SV_ID_0276](scenario_figures/dotplot/PCA_SV_ID_0276.png)
- Haplotype heatmap: [PCA_SV_ID_0276](scenario_figures/heatmap/PCA_SV_ID_0276.png)
- Local-PCA and integrated metrics: [302_region_integrated_metrics.tsv](302_region_integrated_metrics.tsv)

## PCA_SV_ID_0191 — ssa15:51,836,058–53,016,755

This is a secondary candidate. The dotplot contains a long multi-segment pattern with opposite plotted orientations among assemblies and several contig joins. The broad structure is shared, while the detailed arrangement varies. It contains 7 overlapping BayesTyper SVs and 2 direction-discordant comparisons. Relate contains 157 signals and DoLoReS contains 1 signal. The current evidence supports old structural polymorphism with possible lineage-specific rearrangement, but the independent-origin interpretation is weaker than for the three priority candidates above.

- Dotplot: [PCA_SV_ID_0191](scenario_figures/dotplot/PCA_SV_ID_0191.png)
- Haplotype heatmap: [PCA_SV_ID_0191](scenario_figures/heatmap/PCA_SV_ID_0191.png)
- Local-PCA and integrated metrics: [302_region_integrated_metrics.tsv](302_region_integrated_metrics.tsv)

## Required next check

For each priority candidate, compare the actual breakpoint coordinates and orientations in the PAF blocks for each assembly. Independent origin is supported when different lineages have structurally equivalent but breakpoint-distinct configurations, with the same locus-level PCA/haplotype pattern. Shared breakpoints across assemblies support inheritance of an older structural haplotype. The current figures establish candidates for this check; they do not provide the breakpoint test itself.
