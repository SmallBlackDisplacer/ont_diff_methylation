# Identifying regions with different DNA methylation states of homologous chromosomes using Oxford Nanopore reads  

## Description  

Just as single nucleotides, methylation levels may vary on two homologous chromosomes. In nanopore sequencing, electrolytic current signals are sensitive to base modifications, such as 5-methylcytosine (5-mC), so we can determinate methylated positions. The purpose of this study is providing a pipeline to identify areas with varying methylation status on different haplotypes.  

## Methods

We propose a model, in which the number of methylated reads per site is determined by beta binomial distribution. According to null-hypothesis, parameters of beta-binomial distribution are equal for both haplotypes, which alternative hypothesis states that beta-binomial parameters differ. Likelihood is calculated separately for each haplotype and for two haplotypes together. We calculate statistic D = 2·(L1 - L0), which has χ² distribution according to the Wilks’ theorem, and therefore we can calculate p-value that reads on two haplotypes are differentially methylated in a specified genomic window.  

## Dependencies

Nanopolish >= 0.11.3  
Minimap2 >= 2.17-r974-dirty  
Longshot >= 0.4.1  
Snakemake >= 5.19.2  
python >= 3.7   
install the dependencies via  
```pip install -r requirements.txt --user```  

## Analysis workflow examples

```snakemake p_val_result/{sample}.tsv```

## Result example

