# CUT&RUN Pipeline
## Authors: Megan Justice, Alexis Stutzman

## Quick Start:
1. Clone pipeline (Current stable version: v1.0.0)
```
git clone https://github.com/snystrom/cutNrun-pipeline.git --branch v1.0.0 --depth 1 && cd ChIP-pipeline/ && rm -rf .git
```

2. Create `samplesheet.csv` ([see below](#sampleInfo)) with descriptive columns of data.
```
sampleName	fastq	htsfFile
mySample-rep1	htsf_R1.fastq.gz	path/to/htsf_R1.fastq.gz
mySample-rep1	htsf_R2.fastq.gz	path/to/htsf_R2.fastq.gz
```

3. Load python
```
module load python/3.8.8
```

4. Submit the pipeline
```
bash runSnakemake.sh Snakefile.py
```

**Note:** It's a good idea to make sure everything is set up properly before you submit your pipeline. You can do so by first running a "dry run."
```
bash dryrunSnakemake.sh Snakefile.py
```

**Note:** If you're being told there is already an instance of Snakemake happening in your directory.
```
bash dryrunSnakemake.sh Snakefile.py
```