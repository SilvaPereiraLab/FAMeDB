# FAMeDB: A Curated Database for the Analysis of Fungal Aromatic Compound Metabolism 

![banner](https://github.com/SilvaPereiraLab/FAMeDB/blob/main/FAMeDB_Logo.png)

⚠️ FAMeDB will be updated on a regular basis ⚠️

### Introduction



Aromatic compounds represent the second most abundant organic class after carbohydrates, and the study of their microbial metabolism is relevant across several research fields. Numerous peripheral pathways converge on a small number of central intermediates that undergo aromatic ring-opening in the central pathways. This resource is designed to provide a useful and accessible platform for researchers worldwide, even those without specialized expertise in fungal aromatic catabolism, facilitating omics analysis and genomic comparisons.



### Pipeline summary


1. Retrieve single or multiple query proteomes.
2. Run orthology analyis with [proteinortho](https://gitlab.com/paulklemm_PHD/proteinortho).
3. Run Rscript to generate tables and plots.



### Quick Start

1. Download fasta files to directory with FAMeDB.faa file (rename extensiton of fasta to faa if needed)
2. Execute Proteinortho, comparing your .faa file(s) with the FAMeDB: 
```sh
proteinortho directory_fasta_files/*.faa -identity=40 -conn=0.3 -project=Results_FAMeDB
```

Note: For analyses involving 10 or fewer query proteomes, a -conn value of 0.2 may produce improved results.

3. Run the *FAMeDB* script in your terminal or in [RStudio](https://posit.co/download/rstudio-desktop/#download):
```sh
Rscript FAMeDB.R
```

  Note: Results_FAMeDB.proteinortho.tsv file available in [Test with 50 random Ascomycota proteomes](https://github.com/SilvaPereiraLab/FAMeDB/tree/main/Test_output_files)

## Contributions and Support

If you would like to contribute to this database, please see the contributing [guidelines](https://github.com/SilvaPereiraLab/FAMeDB/blob/main/CONTRIBUTING.md).



## Citations


Tiago M. Martins, Rita C. Carmo, Cristina Silva Pereira. 2026. FAMeDB: A Curated Database for the Analysis of Fungal Aromatic Compound Metabolism. Computational and Structural Biotechnology Journal. [doi: 10.34133/csbj.0049](https://doi.org/10.34133/csbj.0049)

Tiago M. Martins, Rita C. Carmo, Cristina Silva Pereira. 2026. FAMeDB: A Curated Database for the Analysis of Fungal Aromatic Compound Metabolism. BioRxiv. [doi: 10.64898/2026.01.19.700319](https://doi.org/10.64898/2026.01.19.700319)



## Acknowledgments

Funding from Fundação para a Ciência e a Tecnologia ([FCT](https://www.fct.pt)) through the project REFILL: Fostering the evolution of fungal biofunnelling of lignin-derived compounds from future biorefineries [doi: 10.54499/2023.15810.PEX](https://doi.org/10.54499/2023.15810.PEX)

TMM is grateful for the working contract (2023.11076.TENURE.076) financed by national funds under the FCT-Tenure Program and RC is grateful to FCT funding for the PhD scholarship 2024.0076.BDANA.

We acknowledge André Cairrão for the design of FAMeDB logo.

