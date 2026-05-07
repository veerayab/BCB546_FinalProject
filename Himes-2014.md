# BCB546_FinalProject: RNA-Seq Transcriptome Profiling Identifies CRISPLD2 as a Glucocorticoid Responsive Gene that Modulates Cytokine Function in Airway Smooth Muscle Cells
Authors: Blanca E. Himes, Xiaofeng Jiang, Peter Wagner, Ruoxi Hu, Qiyu Wang, Barbara Klanderman, Reid M. Whitaker, Qingling Duan, Jessica Lasky-Su, Christina Nikolos, William Jester, Martin Johnson, Reynold A. Panettieri Jr, Kelan G. Tantisira, Scott T. Weiss, Quan Lu

Published: June 13, 2014

Group members: Veeraya Bamrung, Emma Tran

Original paper: https://doi.org/10.1371/journal.pone.0099625

Data GEO: GSE52778

This repository contains R scripts to replicate the RNA-seq analysis from the original paper (Himes et al., 2014) and add new analysis. All scripts are in the 01_scripts folder. The resulting figures and tables are in the 02_output_files folder. The data used are downloaded and placed in the 00_data folder.

## Introduction:
Asthma is a chronic inflammatory lung disease affecting over 25 million Americans and 300 million people worldwide, marked by variable airflow limitation and heightened airway responsiveness. Using RNA-Seq, the authors investigated changes in the transcriptome in four primary human airway smooth muscle (ASM) cell lines that were treated with dexamethasone, an artificial glucocorticoid (GC - common medications used to treat asthma) for cell lines. Results were then validated by functional experiments.

To identify GC-responsive genes in ASM, they performed RNA-seq to profile primary ASM cells from four white male donors treated with dexamethasone. Each sample yielded an average of 58.9 million reads, with about 83.36% aligning to the Homo sapiens (human) genome assembly GRCh37 (hg19) reference genome. Over 98% of mapped bases were mRNA. Transcript and gene expression levels were quantified using hg19 RefSeq annotations. We then used the gene count matrix and reimplemented the analysis with code different from the original paper to replicate Figure 1A and Table 1, since the originally used package (cummeRbund) is outdated and no longer actively maintained. Finally, we performed additional differential gene expression analysis using the DESeq2 package.

## Summary of results:

In conclusion, this study demonstrated how RNA-seq and differential gene expression analysis can be used to investigate the biological effects of dexamethasone treatment. By applying DESeq2 to the RNA-seq data, we identified 316 differentially expressed genes, among which CRISPLD2 emerged as a potential regulator of inflammatory responses. Through this project, we gained practical experience in understanding which types of data are appropriate for DESeq2 analysis and how differential gene expression analysis is performed in practice. At the same time, we learned that reproducing published analyses can be challenging when pipelines and scripts are not thoroughly documented, as reflected by the differences in our results. This project also highlighted the importance of considering version changes in long-term data repositories such as NCBI, as well as the need to carefully evaluate R packages and analysis tools since some may become outdated over time. Overall, this experience provided both biological insights and valuable lessons about the challenges of reproducible bioinformatics research.
