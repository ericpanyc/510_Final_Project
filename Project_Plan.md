# 510_Final_Project

## Study of Alzheimer Disease Correlated Genes

1. The deliverables of this project would be an interactive website accomplished by shiny web app, with which users can view the correlation between AD and several well-studied genes like APP and PS-1. Users can also explore other genes that might be correlated to AD and become potential targets. Finally, users can use the website to see the susceptibility to AD correlated genes' abnormal expression in different parts of the brain.
2. The main dataset used in this project is a table of 50000+ genes' FPKM in 378 samples (both AD and healthy controls, different samples may come from different parts of brain of a same individual). It is a normalized dataset so we could do differential expression analysis on it directly. This dataset, along with other support datasets like donors' clinical information, gene ID conversion table, sample ID conversion table, are provided by Allen Institute for Brain Science, Aging, Dementia and Traumatic Brain Injury Study, and are available at [ADTBIS](http://aging.brain-map.org/overview/home). (Because the main dataset is too large, it can't be uploaded into the github repo)
3. The majority of the project would be differential expression analysis which is done by Limma and Glimma package in R.
4.
    - First Week. Finish dataset processing including removal of lowly-expressed genes, creation of specified sub-dataset for different analysis purpose. Go through vignettes of Limma and Glimma package and finish unit tests.
    - Second Week. Finish differential expression analysis on different datasets created in first week.
    - Third Week. Integrate shiny web app functions to the analysis accomplished in week 2 and release the beta version of the product.
5. The website would let the user choose different genes, then it will display the graphs plotted by R, to help user understand the influence those genes have on AD.
End
