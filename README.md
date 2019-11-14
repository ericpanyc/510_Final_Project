# 510_Final_Project

## Study of Alzheimer Disease Correlated Genes

1. The deliverables of this project would be an interactive website accomplished by shiny web app, with which users can view the correlation between AD and several well-studied genes like APP and PS-1. Users can also explore other genes that might be correlated to AD and become potential targets. Finally, users can use the website to see the susceptibility to AD correlated genes' abnormal expression in different parts of the brain.
2. The data used in this project are 377 RNA-seq BAM files. These data, along with other support datasets like donors' clinical information, gene ID conversion table, sample ID conversion table, are provided by Allen Institute for Brain Science, Aging, Dementia and Traumatic Brain Injury Study, and are available at [ADTBIS](http://aging.brain-map.org/overview/home). (Because the data are too large, it can't be uploaded into the github repo)
3. The majority of the project would be differential expression analysis which is done by Limma and Glimma package in R.
4.
    - First Week. Finish dataset processing including removal of lowly-expressed genes, creation of specified sub-dataset for different analysis purpose. Go through vignettes of Limma and Glimma package and finish unit tests. New milestone: Due to FPKM values are biased and not able to be used in a differential expression analysis, so the first week's plan now is to download all the BAM files available and use figureCounts package to generate counts files.
    - Second Week. Finish differential expression analysis on different datasets created in first week.
    - Third Week. Integrate shiny web app functions to the analysis accomplished in week 2 and release the beta version of the product.
5. The website would let the user choose different genes, then it will display the graphs plotted by R, to help user understand the influence those genes have on AD.
End
