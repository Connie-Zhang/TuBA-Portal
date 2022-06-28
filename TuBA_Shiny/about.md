## About

This app demonstrates using R shiny to visualize results of the TuBA algorithm. It offers graphs on the survival curve and copy number of biclusters for 23 cancer types. More importantly, this application allows the filtering of biclusters based on genes, biological pathways, significance of survival, and copy number enrichment. 

In the “Biological Pathways” page, the "Select the bicluster of interest" input is adjusted based on the specific gene and biological pathway of interest. Users can also choose to filter by significance of biological pathway and survival analysis. In the Gene Information tab, a table of gene information is also provided, including gene names and chromosome. A vector of all the genes in the bicluster is also provided. 

In the “Copy Number” page, the "Select the bicluster of interest" input is adjusted based on the chromosome of interest. It provides a filtered list of biclusters where more than 80% of its genes belong in the selected chromosome of interest. The gene information table and vector are provided as well. 

If wish to compare the survival curve and copy number enrichment of a bicluster, or compare different cancer types, another website tab can be opened for comparison. 

## Credits

This project was developed by [Connie Zhang] (https://connie-zhang.github.io) while she was an intern at the [Khiabanian Lab] (http://www.khiabanian-lab.org), as a part of the Rutgers University DIMACS REU program. Amartya Singh, is a postdoc in our lab, developed the TuBA algorithm and provided useful input in how to employ the algorithm results for data analysis. 

The plotting design is developed using the ggplot2 package. The TCGA expression data were downloaded from the [USCS Xena portal] (https://xenabrowser.net/datapages/). 

Authors:
- Connie Zhang
- Amartya Singh
- Hossein Khiabanian