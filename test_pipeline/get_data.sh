# Get test data!
wget http://s3.climb.ac.uk/nanopolish_tutorial/methylation_example.tar.gz
tar -xvf methylation_example.tar.gz

# And run snakemake!
snakemake add_ht/albacore_output.tsv
