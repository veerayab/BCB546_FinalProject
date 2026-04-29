# BCB546_FinalProject: RNA-Seq Differential Expression
Group members: Veeraya Bamrung, Emma Tran

Original paper: https://doi.org/10.1371/journal.pone.0099625

Data GEO: GSE52778

This repository contains R scripts to replicate the RNA-seq analysis from the original paper (Himes et al., 2014) and add new analysis. All scripts are in the 01_scripts folder. The resulting figures and tables are in the 02_output_files folder. The data used are downloaded and placed in the 00_data folder.

Asthma is a chronic inflammatory lung disease affecting over 25 million Americans and 300 million people worldwide, marked by variable airflow limitation and heightened airway responsiveness. Using RNA-Seq, the authors investigated changes in the transcriptome in four primary human airway smooth muscle (ASM) cell lines that were treated with dexamethasone, an artificial glucocorticoid (GC - common medications used to treat asthma) for cell lines. Results were then validated by functional experiments.

To identify GC-responsive genes in ASM, they performed RNA-seq to profile primary ASM cells from four white male donors treated with dexamethasone. Each sample yielded an average of 58.9 million reads, with about 83.36% aligning to the Homo sapiens (human) genome assembly GRCh37 (hg19) reference genome. Over 98% of mapped bases were mRNA. Transcript and gene expression levels were quantified using hg19 RefSeq annotations. We then used the gene count matrix and reimplemented the analysis with code different from the original paper to replicate Figure 1A and Table 1, since the originally used package (cummeRbund) is outdated and no longer actively maintained. Finally, we performed additional differential gene expression analysis using the DESeq2 package.
