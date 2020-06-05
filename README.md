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

+------+-------+-------+-------------+-------------+---------+-------+-------+---------+-------+-------+----------+-------+-------+--------+
|chrom | start | end   | a1          | b1	       | l1      | a2    | b2    | l2      | a0    | b0    | l0       | D     | n_pos | p_val  |
+------+-------+-------+-------------+-------------+---------+-------+-------+---------+-------+-------+----------+-------+-------+--------+
|chr9  | 35148 | 37148 | 8337959.275 | 8337809.099 | -5.628  | 0.709 | 1.024 | -8.378  | 1.417 | 1.771 | -15.751  | 3.490 | 3     | 0.1747 |
+------+-------+-------+-------------+-------------+---------+-------+-------+---------+-------+-------+----------+-------+-------+--------+
|chr9  | 37148 | 39148 | 1.144       | 2.033       | -23.943 | 1.112 | 2.863 | -24.792 | 1.090 | 2.330 | -49.016  | 0.561 | 10    | 0.7556 |
+------+-------+-------+-------------+-------------+---------+-------+-------+---------+-------+-------+----------+-------+-------+--------+
|chr9  | 39148 | 41148 | 1.053       | 2.276       | -50.388 | 1.053 | 2.923 | -53.596 | 1.039 | 2.544 | -104.246 | 0.523 | 21    | 0.7699 |
+------+-------+-------+-------------+-------------+---------+-------+-------+---------+-------+-------+----------+-------+-------+--------+
|chr9  | 41148 | 43148 | 0.743       | 1.074       | -32.217 | 0.825 | 1.513 | -32.759 | 0.771 | 1.250 | -65.146  | 0.341 | 12    | 0.8433 |
+------+-------+-------+-------------+-------------+---------+-------+-------+---------+-------+-------+----------+-------+-------+--------+
|chr9  | 43148 | 45148 | 1.110       | 1.262       | -23.609 | 0.886 | 2.037 | -19.529 | 0.926 | 1.457 | -43.857  | 1.437 | 8     | 0.4874 |
+------+-------+-------+-------------+-------------+---------+-------+-------+---------+-------+-------+----------+-------+-------+--------+
