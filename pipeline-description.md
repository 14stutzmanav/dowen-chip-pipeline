# Dowen Lab ChIP Pipeline Description
This document is for people just beginning their journey through analyzing data. First and foremost, welcome! Computational biology is a fun field and I'm excited that you've decided to give it a shot. In this document, we're going to walk through the basics of our ChIP pipeline and point out some areas you want to pay attention to. Note that this is a description of the pipeline steps, not a description of how snakemake works generally. If you're curious about that, check out this resource: https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html


## Part 1: Housekeeping rules
Okay, let's begin our walk through of the pipeline. The first rule is called `rule all`. This rule is a set up for the rest of the pipeline and it gives Snakemake the ability to pick up where it last left off if your pipeline geets interrupted during a run. Basically, it just "expands" the subset files produced by the pipeline that are needed for subsequent rules. Unless you're editing the pipeline structure, you most likely don't need to edit this step. For ease of reference, here's the complete `rule all`:
```
rule all:
	input:
		expand('Fastq/{fastq}', fastq = list(sampleDF.fastq)),
		expand('Fastq/{sample}_R{Num}.fastq.gz', sample = sampleList, Num = ['1', '2']),
		expand('FastQC/{sample}_R1_fastqc.html', sample = sampleList),
		expand('Sam/{sample}.sam', sample = sampleList),
		expand('Bam/{sample}-fixmate.bam', sample = sampleList),
		expand('Bam/{sample}-fixmate_sort.bam', sample = sampleList),
                expand('Bam/{sample}-notfixed-sort.bam', sample = sampleList),
		expand('Bam/{sample}-notfixed-sortedbyname.bam', sample = sampleList),
                expand('Bam/{sample}-markdup.bam', sample = sampleList),
		expand('Bam/{sample}-markdup-sort.bam', sample = sampleList),
                expand('Bam/{sample}-markdup-sort.bam.bai', sample = sampleList),
		expand('Counts/{sample}-mousealign.txt', sample = sampleList),
                expand('Bam/{sample}-mousereads.bam', sample = sampleList),
                expand('Bam/{sample}-mousereads-fix.bam', sample = sampleList),
		expand('Bam/{sample}-mousereads-fix-sort.bam', sample = sampleList),
		expand('Bed/{sample}-fixall.bed', sample = sampleList),
		expand('Counts/{sample}-NormFactor.txt', sample = sampleList),
		expand('Bedgraph/{sample}-sort.bedgraph', sample = sampleList),
        	expand('Bigwig/{sample}.bw', sample = sampleList)
```

Next, we have `rule copyFiles`. This rule is going to look at your sampleSheet and copy fastq files from `HTSF` to a new folder, called `Fastq`, in your current directory (I call this your working directory, since it's where you're "working" from). When it's running, there should be a little terminal message that says `Copying files to Fastq directory with corrected file names.` Here's what that rule looks like:
```
rule copyFiles:
	input:
		lambda x: list(sampleDF.htsfFile)
	output:
		expand('Fastq/{fastq}', fastq = list(sampleDF.fastq))
	message: "Copying files to Fastq directory with corrected file names."
        run:
		for htsf in list(sampleDF.htsfFile):
			outFileFilt = sampleDF [ sampleDF.htsfFile == htsf ] 
			outFileBase = list(outFileFilt.fastq)[0].split('.')[0]
			outFile = 'Fastq/{fastq}.fastq.gz'.format(fastq = outFileBase)
			shutil.copyfile(htsf, outFile)
			print('copied file')

```

The last "housekeeping step" for the pipeline is `combinesampleReps`. 
