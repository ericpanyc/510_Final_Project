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


## Milestone Update:

1. Milestone 1 (11/14/18): Found that FPKM value is not suitable for differential expression analysis, need to use bam data to generate raw counts table. Following step: Download all bam files and learn how to use featureCounts function.
2. Milestone 2 (11/21/18): Could get a raw counts table of gene counts matrix from all samples, but the DGEList I got is not identical to the sample on the vignette. Following step: Need to figure out how to follow the vignette to proceed to the next step.
3. Milestone 3 (11/27/18): Could get output images following the limma vignette. Following step: Need to set up an appropriate design matrix which can best reflect the genes' correlation with AD.
4. Final Update. The application for PCA was published. This application uses RNA lcpm data as input, to do a PCA on AD (whether diagnosed with AD) or sex factor. But it does not function consistently (could only work when select AD as input value) and the PCA result is not good.

## Project Conclusion:
###  As still an untreatable disease, Alzheimer Disease is extremely complicated at RNA level. Using basic analysis method like PCA could not get satisfactory results. So, more advanced methods in the research of AD are in urgent need.  

- Published R Shiny app [ADPlot](https://ericpanyc.shinyapps.io/ADPlot/)
- Published R notebook [Project Notebook](http://rpubs.com/ericpanyc/557398)
- Large Files could be downloaded here [Large Files](https://drive.google.com/drive/folders/1crh2GuOxyJKxGpncZWn_TXzn9oMQtceJ)
