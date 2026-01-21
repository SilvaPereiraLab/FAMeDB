# FAMeDB: A Curated Database for the Analysis of Fungal Aromatic Compound Metabolism 

![banner](https://github.com/SilvaPereiraLab/FAMeDB/blob/main/FAMeDB_Logo.png)


### Introduction



Aromatic compounds represent the second most abundant organic class after carbohydrates, and the study of their microbial metabolism is relevant across several research fields. Numerous peripheral pathways converge on a small number of central intermediates that undergo aromatic ring-opening in the central pathways. This resource is designed to provide a useful and accessible platform for researchers worldwide, even those without specialized expertise in fungal aromatic catabolism, facilitating omics analysis and genomic comparisons.



### Pipeline summary



1. Run orthology analyis with [proteinortho](https://gitlab.com/paulklemm_PHD/proteinortho).
2. Run Rscript to generate tables and plots.



### Quick Start

1. Download fasta files to directory with FAMeDB.faa file (rename extensiton of fasta to faa if needed)
2. Execute Proteinortho, comparing your .faa file(s) with the FAMeDB: 
```sh
proteinortho directory_fasta_files/*.faa -identity=40 -conn=0.3 -project=Results_FAMeDB
```
3. Run the *FAMeDB* script in your terminal or in [RStudio](https://posit.co/download/rstudio-desktop/#download):
```sh
Rscript FAMeDB.R
```


## Contributions and Support

If you would like to contribute to this database, please see the contributing [guidelines](link).



## Citations


(Under peer-review submission process)

Tiago M. Martins, Rita C. Carmo, Cristina Silva Pereira. 2026. FAMeDB: A Curated Database for the Analysis of Fungal Aromatic Compound Metabolism. BioRxiv. [doi: 10.64898/2026.01.19.700319](https://doi.org/10.64898/2026.01.19.700319)



## Acknowledgments

Funding from Fundação para a Ciência e a Tecnologia ([FCT](https://www.fct.pt)) through the project REFILL: Fostering the evolution of fungal biofunnelling of lignin-derived compounds from future biorefineries [doi: 10.54499/2023.15810.PEX](https://doi.org/10.54499/2023.15810.PEX)

TMM is grateful for the working contract (2023.11076.TENURE.076) financed by national funds under the FCT-Tenure Program and RC is grateful to FCT funding for the PhD scholarship 2024.0076.BDANA.

We acknowledge André Cairrão for the design of FAMeDB logo.

