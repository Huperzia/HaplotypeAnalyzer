# HaplotypeAnalyzer
A shiny app for haplotype analysis in a diverse population of sweet corn

This app allows you to input a location on the Ia453-sh2 sweet corn genome and will perform clustering and plotting of the haplotypes at that locus. Additionally, accession IDs can be input to color populations and phenotypes can be input to plot alongside the haplotypes to visualize and mine associations. 

Below are screenshots of different loci. Gray is ref allele, green is alternate allele, orange is heterozygote, and purple is missing (with high enough depth missing is likely an SV - particularly in the se1 allele)

Haplotype analysis of shrunken2 in the sweet corn diversity population: 

![su1_haplotype](./figs/sh2_haplotype.png)

Haplotype analysis of sugary1 in the sweet corn diversity population:

![su1_haplotype](./figs/su1_haplotype.png)

Haplotype analysis of sugary enhancer1 in the sweet cap population:

![su1_haplotype](./figs/se1_haplotype.png)

Additional features include on-the-fly pca and dapc analyses for any given locus. Here is an example of a sh2 locus with mutations highlighted.
![hap_pca](./figs/hap_pca.png)
