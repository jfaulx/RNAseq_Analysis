# Welcome! 
## This is an RShiny app that I created to help analyze and visualize trends in RNA sequencing data.

Features include metadata summary, count data filtering and visualization, differential expression via DESeq2, and individual gene expression analysis
 
You can find the app hosted here: https://jacksonfaulx.shinyapps.io/rnaseqapp/ 

The app is still being actively worked on and I am looking to add new features regularly. It currently allows you to upload a RNAseq raw counts file as well as a meta data file, formatted as seen below. Both files can be .csv or .tsv


Counts file: Needs to have gene IDs as row names and samples as columns, counts data should be non-normalized and non-filtered

![image](https://github.com/jfaulx/RNAseq_Analysis/assets/143756015/01ffb169-b3eb-43b3-93fa-6e0eb6043864)


Meta data: Should have columns containing sample information, first column must be called 'sample' and contain all sample names

![image](https://github.com/jfaulx/RNAseq_Analysis/assets/143756015/736e8579-8148-4e07-a691-6aa108e5b155)

Example files can also be found in this repo, they are derived from a study comparing cardiac myocytes in neonatal and adult mice

(O'Meara CC, Wamstad JA, Gladstone RA, Fomovsky GM et al. Transcriptional reversion of cardiac myocyte fate during mammalian cardiac regeneration. Circ Res 2015 Feb 27;116(5):804-15. PMID: 25477501)

## Thank you! I hope you enjoy
