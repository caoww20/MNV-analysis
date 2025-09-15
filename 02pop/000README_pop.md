## MNV overlap in populations
## Impact of MNVs on population stratification (using old data, as new data includes MHC, black, and X, but does not affect this)

## Plotted images
1. Venn diagram of MNVs in populations
2. PCA plots of SNVs and MNVs
    SNV and MNV PCA clustering effects are similar, incidentally check the heritability of variants
3. FST Manhattan plots, correlation plots, population heatmaps, population heatmap differences for SNPs and MNVs
    Since FST uses genotype frequencies to measure genetic differentiation, here we use FST instead of genotype frequencies
    Calculate FST for SNPs and MNVs by interval, draw Manhattan plots, and calculate their correlation
    Draw heatmaps between populations, including SNV, MNV, and their differences: compared to SNV, MNV has higher genetic differentiation ratios in which populations, which is to prove that although MNV is very similar to SNV, there are still differences in populations
4. Bar chart of the number of MNVs with FST > 0.15 in different populations
    Calculate the number of MNVs where FST is greater than 0.15 compared to any one population, indicating whether a certain population has more such differences, or all have them
5. Enrichment analysis plot
    Taking AFR as an example, map it to genes, see the results of enrichment analysis【See if it is related to weather, African diseases, etc.】
    # Later discovered that using all data with FST >= 0.15 sites for enrichment analysis, genes with at least 3 significant FST variants are included in enrichment analysis
    Mainly concentrated in cardiomyopathy, understand the epidemiology of populations as annotation
    https://pubmed.ncbi.nlm.nih.gov/16330699/
    https://pubmed.ncbi.nlm.nih.gov/23536465/
    https://www.nature.com/articles/s41569-020-0428-2
    (Arrhythmogenic right ventricular cardiomyopathy[Title/Abstract]) AND (Epidemiology[Title/Abstract])
    (cardiomyopathy[Title/Abstract]) AND (Epidemiology[Title/Abstract])
    # View which genes have the highest mutations

## Whether African-specific or high-frequency MNVs are related to UV intensity!!!
## Comparison between African populations and European or East Asian populations, see if some mutations are concentrated elsewhere
# Population-specific
    # Perform enrichment analysis on specific mutations
# Population-shared
    # Use FST, genotype clustering, etc. for population clustering analysis
    # Draw outlier loci based on MNV frequencies in different populations (using FST instead)
# Specific and outlier mean evolution
# Draw phylogenetic trees (can be omitted)
## Population analysis (temporarily abandoned)
    ## Two strategies
    ## First strategy according to mutation pattern frequency
        ## First splice 1000G mnv and mnv_multi from luohh's folder

        ## Use 01filter1000G.py to filter points with more than 5 joints, seems not done?
        ## Use MNVAnno to annotate mnv.txt to get mnv.genotype with 0,1,2 genotypes
        ## Use 02group.py to get the frequency of each type in each individual
        ## Use 00plot.R for principal component analysis
    ## Second strategy calculate according to FST
        ## Use 01getFST.py to get FST merged results with ACB as core
        ## Use 00plot.R for principal component analysis


