# Identifying regions with different DNA methylation states of homologous chromosomes using Oxford Nanopore reads  

## Description  

Just as single nucleotides, methylation levels may vary on two homologous chromosomes. At the same time, DNA methylation is now recognised as important epigenetic mechanism that plays a crucial role in cellular regulatory systems. Recently DNA methylation variation has been proposed to be a potential contributing factor to cancer risk. In nanopore sequencing, electrolytic current signals are sensitive to base modifications, such as 5-methylcytosine (5-mC), so we can determinate methylated positions.   
The purpose of this study is providing a pipeline to identify areas with varying methylation status on different haplotypes. Methylated CpG sites were determined using Nanopolish, and each read was assigned to its own haplotype using long-read variant caller LongShot. Then, we can estimate the number of methylated and non-methylated reads for each site and each haplotype.  
We propose a model, in which the number of methylated reads per site is determined by beta binomial distribution. According to null-hypothesis, parameters of beta-binomial distribution are equal for both haplotypes, which alternative hypothesis states that beta-binomial parameters differ. Likelihood is calculated separately for each haplotype and for two haplotypes together. We calculate statistic D = 2·(L1 - L0), which has χ² distribution according to the Wilks’ theorem, and therefore we can calculate p-value that reads on two haplotypes are differentially methylated in a specified genomic window.  
Regions that show multiple sites with significant methylation differences can be considered differentially methylated. The method was applied to whole genome data for the HG001/NA12878 individual from the Nanopore WGA Consortium.  

## Pipeline



## Dependencies

python >= 3.7  
install the dependencies via  
```pip install -r scripts/requirements.txt --user```  

## Analysis workflow examples
