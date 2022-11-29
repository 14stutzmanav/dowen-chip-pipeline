## ChIP pipeline was written by Megan Justice and adapted to Snakemake by Alexis Stutzman
## Last edited Nov 2022.

import pandas as pd
import shutil
import os
import re
import glob
import getpass
import sys
import argparse

#################################################################################
# Startup Instructions (check out the READ.me for more details!):
# 1. Before you begin, make sure you update 'example-samplesheet.csv' in "sampleSheetPath = str('example-samplesheet.csv')" to point to your .csv
#
# 2. If needed, edit Genome_DIR, chromSizes, and REPEATPATH to point to appropriate files.
#    -> Note that you only need to do this if you're changing the assembly.
#
# 3. Run the pipeline. First you need to load python (module load python), then start the pipeline (bash runSnakemake.sh Snakefile.py)
#    -> It's good practice to make sure everything is set up properly by first doing a "dry run." (bash dryrunSnakemake.sh Snakefile.py)
#    -> If the computer thinks you're running an instance of snakemake already (and you've checked that no jobs from the pipeline are 
#          currently running or pending in SLURM), you need to "unlock" them  directory (bash unlockSnakemake.sh Snakefile.py)
#
# Troubleshooting:
# If a step needs more time or memory, edit slurmConfig.json
#################################################################################
sampleSheetPath = str('example-samplesheet.csv')

# Alex's TODO: Move to config file (see "TODO below")
### General Dowen Lab setup:
Genome_DIR = str('/proj/dowenlab/projects/ChIP-seq/MergedGenomeforNormalization/MergedGenome')
chromSizes = str('/proj/dowenlab/projects/ChIP-seq/mm10.chrom.sizes')
REPEATPATH = str('/proj/dowenlab/projects/ChIP-seq/RepeatSequence.bed')

### Audra's setup:
# Genome_DIR=/proj/dowenlab/users/megan/MergedGenomeforNormalization/MergedGenome
# chromSizes='/proj/dowenlab/users/megan/Scripts/SpikeInNormalization/mm10.chrom.sizes'
# DIR=/proj/dowenlab/users/audra/siGlo-ChIP/221012_UNC51-VH00562_132_AAATJ35HV/Peaks
# REPEATPATH = RepeatSequence.bed

#################################################################################
# Pipeline setup (don't change this part unless you REALLY have to)
#################################################################################

configfile: 'config.json'
modules = config['module']

sampleDF = pd.read_csv(sampleSheetPath, comment = '#')
sampleList = list(set(sampleDF.sampleName))

def getsampleReplicates(wildcards):
        readNumRegex = '_R{}'.format(wildcards.Num)
        sampleFilter = sampleDF[sampleDF.sampleName == wildcards.sample]
        fastqList = list(sampleFilter.fastq)
        fastqList = [ fastq for fastq in fastqList if re.search(readNumRegex, fastq) ]            
        fastqList = [ 'Fastq/{}'.format(fastq) for fastq in fastqList ]

        return(fastqList)

#constrain the sample wildcard to characters, numbers, or underscore
#wildcard_constraints:
#    sample = '\w+'

# TODO: Revise ref and spike genome to make "MergedGenomeforNormalization" each time the pipeline is run. This gives users flexibility to change ref and spike genomes easily.
# REFGENOME = config['refGenome']
# SPIKEGENOME = config['spikeGenome']
# REFGENOMEPATH = config['genome'][REFGENOME]['bowtie']
# chromSize_Path  = config['genome'][REFGENOME]['chrSize']
# genomeSize = config['genome'][REFGENOME]['genomeSize']
# readLen = config['readLen']


#################################################################################
# Begin pipeline
#################################################################################

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

rule combinesampleReps:
	input:
		getsampleReplicates
	output:
                'Fastq/{sample}_R{Num,[12]}.fastq.gz'
	shell:
		"""
		cat {input} > {output}
		"""

rule fastQC:
	input:
		'Fastq/{sample}_R1.fastq.gz',
	output:
		'FastQC/{sample}_R1_fastqc.html'
	envmodules:
		modules['fastqcVer']
	benchmark:
		"benchmarks/{sample}.fastQC.benchmark.txt"
	shell:
		"""
		fastqc -o ./FastQC/ -f fastq {input}
		"""

rule align:
	input:
		r1 = 'Fastq/{sample}_R1.fastq.gz',
		r2 = 'Fastq/{sample}_R2.fastq.gz'
	output:
		sam = 'Sam/{sample}.sam',
		logInfo = 'Logs/{sample}.log'
	message: "Aligning fastq files. This could take a while..."
	params:
		genome = Genome_DIR
	envmodules:
		modules['bowtieVer']
	shell:
		"""
		(bowtie -v 2 -p 12 -S --seed 123 -m 1 {params.genome} -1 {input.r1} -2 {input.r2} > {output.sam}) 2> {output.logInfo}
		"""

rule convertToBamAndFix:
	input:
		'Sam/{sample}.sam'
	output:
		unsorted = 'Bam/{sample}.bam',
		notfixed = 'Bam/{sample}-notfixed-sort.bam',
		namesort = 'Bam/{sample}-notfixed-sortedbyname.bam',
		fixmate = 'Bam/{sample}-fixmate.bam',
		fixsort = 'Bam/{sample}-fixmate_sort.bam',
		dupsMarked = 'Bam/{sample}-markdup.bam',
		dupsSort = 'Bam/{sample}-markdup-sort.bam',
		index = 'Bam/{sample}-markdup-sort.bam.bai',
		mousereads = 'Bam/{sample}-mousereads.bam',
		mousefix = 'Bam/{sample}-mousereads-fix.bam',
		mousesort = 'Bam/{sample}-mousereads-fix-sort.bam'
	threads: 12
	message: "Congratulations! Your files aligned. Now converting .sam to .bam and sorting the .bam file."
	envmodules:
		modules['samtoolsVer']
	shell:
		"""
		samtools view -b -S -F 4 -@ {threads} -o {output.unsorted} {input}
		samtools sort -@ {threads} -o {output.notfixed} {output.unsorted}
		samtools sort -n -@ {threads} {output.notfixed} -o {output.namesort}
		samtools fixmate -m -@ {threads} {output.namesort} {output.fixmate}
		samtools sort -@ {threads} {output.fixmate} -o {output.fixsort}
		samtools markdup -r -s -@ {threads} {output.fixsort} {output.dupsMarked}
		samtools sort -@ {threads} -o {output.dupsSort} {output.dupsMarked}
		samtools index {output.dupsSort} {output.index}
		samtools view -h -@ {threads} {output.dupsSort} | grep Mchr > {output.mousereads}
		sed 's/Mchr/chr/g' {output.mousereads} > {output.mousefix}
		samtools sort {output.mousefix} > {output.mousesort}
		"""

rule countMouse:
	input:
		'Bam/{sample}-markdup-sort.bam'
	output:
		mousealign = 'Counts/{sample}-mousealign.txt',
		mousecount = 'Counts/{sample}-mousecount.txt'
	envmodules:
		modules['samtoolsVer']
	shell:
		"""
		samtools idxstats {input} | grep Mchr > {output.mousealign}
		awk '{{s+=$3}}END{{print s}}' {output.mousealign} > {output.mousecount}
		"""

rule countHuman:
	input:
		'Bam/{sample}-markdup-sort.bam'
	output:
		humalign = 'Counts/{sample}-humanalign.txt',
		humcount = 'Counts/{sample}-humancount.txt',
		normFactor = 'Counts/{sample}-NormFactor.txt'
	envmodules:
		modules['samtoolsVer']
	shell:
		"""
		samtools idxstats {input} | grep -v Mchr > {output.humalign}
		awk '{{s+=$3}}END{{print s}}' {output.humalign} > {output.humcount}
		normfactor=$( awk '{{h+=$1}}END{{if(h<50000){{print 1}} else {{print (1 / (h/1000000))}}}}' {output.humcount} )
		echo $normfactor > {output.normFactor}
		"""

rule bamToBed:
	input:
		'Bam/{sample}-mousereads-fix-sort.bam'
	output:
		bed = 'Bed/{sample}.bed',
		extended = 'Bed/{sample}-sort-extend.bed',
		fixstart = 'Bed/{sample}-fixstart.bed',
		fixall = 'Bed/{sample}-fixall.bed'
	envmodules:
		modules['bedtoolsVer']
	shell:
		"""
		bedtools bamtobed -i {input} > {output.bed}
		awk '{{
			if ($6 == "-")
				print $1, $2-200, $3, $4, $5, $6;
			else
				print $1, $2, $3+200, $4, $5, $6;
			}}' {output.bed} > {output.extended}

		awk '{{
			if ($2 < 0)
				print $1, 0, $3, $4, $5, $6;
			else
				print $1, $2, $3, $4, $5, $6;
			}}' {output.extended} > {output.fixstart}

		awk '{{OFS="\t"}}{{
			if (($1 == "chr1") && ($3 > 195471971))
				print $1, $2, 195471971, $4, $5, $6;
			else if (($1 == "chr2") && ($3 > 182113224))
				print $1, $2, 182113224, $4, $5, $6;
			else if (($1 == "chr3") && ($3 > 160039680))
				print $1, $2, 160039680, $4, $5, $6;
			else if (($1 == "chr4") && ($3 > 156508116))
				print $1, $2, 156508116, $4, $5, $6;
			else if (($1 == "chr5") && ($3 > 151834684))
				print $1, $2, 151834684, $4, $5, $6;
			else if (($1 == "chr6") && ($3 > 149736546))
				print $1, $2, 149736546, $4, $5, $6;
			else if (($1 == "chr7") && ($3 > 145441459))
				print $1, $2, 145441459, $4, $5, $6;
			else if (($1 == "chr8") && ($3 > 129401213))
				print $1, $2, 129401213, $4, $5, $6;
			else if (($1 == "chr9") && ($3 > 124595110))
				print $1, $2, 124595110, $4, $5, $6;
			else if (($1 == "chr10") && ($3 > 130694993))
				print $1, $2, 130694993, $4, $5, $6;
			else if (($1 == "chr11") && ($3 > 122082543))
				print $1, $2, 122082543, $4, $5, $6;
			else if (($1 == "chr12") && ($3 > 120129022))
				print $1, $2, 120129022, $4, $5, $6;
			else if (($1 == "chr13") && ($3 > 120421639))
				print $1, $2, 120421639, $4, $5, $6;
			else if (($1 == "chr14") && ($3 > 124902244))
				print $1, $2, 124902244, $4, $5, $6;
			else if (($1 == "chr15") && ($3 > 104043685))
				print $1, $2, 104043685, $4, $5, $6;
			else if (($1 == "chr16") && ($3 > 98207768))
				print $1, $2, 98207768, $4, $5, $6;
			else if (($1 == "chr17") && ($3 > 94987271))
				print $1, $2, 94987271, $4, $5, $6;
			else if (($1 == "chr18") && ($3 > 90702639))
				print $1, $2, 90702639, $4, $5, $6;
			else if (($1 == "chr19") && ($3 > 61431566))
				print $1, $2, 61431566, $4, $5, $6;
			else if (($1 == "chrX") && ($3 > 171031299))
				print $1, $2, 171031299, $4, $5, $6;
			else if (($1 == "chrY") && ($3 > 91744698))
				print $1, $2, 91744698, $4, $5, $6;
			else if (($1 == "chrM") && ($3 > 16299))
				print $1, $2, 16299, $4, $5, $6;
			else if (($1 == "chr5_JH584299_random") && ($3 > 953012))
				print $1, $2, 953012, $4, $5, $6;
			else if (($1 == "chrX_GL456233_random") && ($3 > 336933))
				print $1, $2, 336933, $4, $5, $6;
			else if (($1 == "chrY_JH584301_random") && ($3 > 259875))
				print $1, $2, 259875, $4, $5, $6;
			else if (($1 == "chr1_GL456211_random") && ($3 > 241735))
				print $1, $2, 241735, $4, $5, $6;
			else if (($1 == "chr4_GL456350_random") && ($3 > 227966))
				print $1, $2, 227966, $4, $5, $6;
			else if (($1 == "chr4_JH584293_random") && ($3 > 207968))
				print $1, $2, 207968, $4, $5, $6;
			else if (($1 == "chr1_GL456221_random") && ($3 > 206961))
				print $1, $2, 206961, $4, $5, $6;
			else if (($1 == "chr5_JH584297_random") && ($3 > 205776))
				print $1, $2, 205776, $4, $5, $6;
			else if (($1 == "chr5_JH584296_random") && ($3 > 199368))
				print $1, $2, 199368, $4, $5, $6;
			else if (($1 == "chr5_GL456354_random") && ($3 > 195993))
				print $1, $2, 195993, $4, $5, $6;
			else if (($1 == "chr4_JH584294_random") && ($3 > 191905))
				print $1, $2, 191905, $4, $5, $6;
			else if (($1 == "chr5_JH584298_random") && ($3 > 184189))
				print $1, $2, 184189, $4, $5, $6;
			else if (($1 == "chrY_JH584300_random") && ($3 > 182347))
				print $1, $2, 182347, $4, $5, $6;
			else if (($1 == "chr7_GL456219_random") && ($3 > 175968))
				print $1, $2, 175968, $4, $5, $6;
			else if (($1 == "chr1_GL456210_random") && ($3 > 169725))
				print $1, $2, 169725, $4, $5, $6;
			else if (($1 == "chrY_JH584303_random") && ($3 > 158099))
				print $1, $2, 158099, $4, $5, $6;
			else if (($1 == "chrY_JH584302_random") && ($3 > 155838))
				print $1, $2, 155838, $4, $5, $6;
			else if (($1 == "chr1_GL456212_random") && ($3 > 153618))
				print $1, $2, 153618, $4, $5, $6;
			else if (($1 == "chrUn_JH584304") && ($3 > 114452))
				print $1, $2, 114452, $4, $5, $6;
			else if (($1 == "chrUn_GL456379") && ($3 > 72385))
				print $1, $2, 72385, $4, $5, $6;
			else if (($1 == "chr4_GL456216_random") && ($3 > 66673))
				print $1, $2, 66673, $4, $5, $6;
			else if (($1 == "chrUn_GL456393") && ($3 > 55711))
				print $1, $2, 55711, $4, $5, $6;
			else if (($1 == "chrUn_GL456366") && ($3 > 47073))
				print $1, $2, 47073, $4, $5, $6;
			else if (($1 == "chrUn_GL456367") && ($3 > 42057))
				print $1, $2, 42057, $4, $5, $6;
			else if (($1 == "chrUn_GL456239") && ($3 > 40056))
				print $1, $2, 40056, $4, $5, $6;
			else if (($1 == "chr1_GL456213_random") && ($3 > 39340))
				print $1, $2, 39340, $4, $5, $6;
			else if (($1 == "chrUn_GL456383") && ($3 > 38659))
				print $1, $2, 38659, $4, $5, $6;
			else if (($1 == "chrUn_GL456385") && ($3 > 35240))
				print $1, $2, 35240, $4, $5, $6;
			else if (($1 == "chrUn_GL456360") && ($3 > 31704))
				print $1, $2, 31704, $4, $5, $6;
			else if (($1 == "chrUn_GL456378") && ($3 > 31602))
				print $1, $2, 31602, $4, $5, $6;
			else if (($1 == "chrUn_GL456389") && ($3 > 28772))
				print $1, $2, 28772, $4, $5, $6;
			else if (($1 == "chrUn_GL456372") && ($3 > 28664))
				print $1, $2, 28664, $4, $5, $6;
			else if (($1 == "chrUn_GL456370") && ($3 > 26764))
				print $1, $2, 26764, $4, $5, $6;
			else if (($1 == "chrUn_GL456381") && ($3 > 25871))
				print $1, $2, 25871, $4, $5, $6;
			else if (($1 == "chrUn_GL456387") && ($3 > 24685))
				print $1, $2, 24685, $4, $5, $6;
			else if (($1 == "chrUn_GL456390") && $3 > 24668))
				print $1, $2, 24668, $4, $5, $6;
			else if (($1 == "chrUn_GL456394") && ($3 > 24323))
				print $1, $2, 24323, $4, $5, $6;
			else if (($1 == "chrUn_GL456392") && ($3 > 23629))
				print $1, $2, 23629, $4, $5, $6;
			else if (($1 == "chrUn_GL456382") && ($3 > 23158))
				print $1, $2, 23158, $4, $5, $6;
			else if (($1 == "chrUn_GL456359") && ($3 > 22974))
				print $1, $2, 22974, $4, $5, $6;
			else if (($1 == "chrUn_GL456396") && ($3 > 21240))
				print $1, $2, 21240, $4, $5, $6;
			else if (($1 == "chrUn_GL456368") && ($3 > 20208))
				print $1, $2, 20208, $4, $5, $6;
			else if (($1 == "chr4_JH584292_random") && ($3 > 14945))
				print $1, $2, 14945, $4, $5, $6;
			else if (($1 == "chr4_JH584295_random") && ($3 > 1976))
				print $1, $2, 1976, $4, $5, $6;
			else
				print $1, $2, $3, $4, $5, $6;
		}}' {output.fixstart} > {output.fixall}
		"""

rule callPeaks:
	input:
		'Bed/{sample}-fixall.bed'
	output:
		summits = 'Bed/{sample}_summits.bed',
		exp50 = 'Bed/{sample}-summits-exp50.bed',
		noreps = 'Bed/{sample}-summits-exp50-noreps.bed',
		shrink = 'Bed/{sample}-summits-exp50-noreps-shrink.bed'
	message: "The pipeline is now calling peaks and fixing summit files ('Hi, is this Peaks? Yeah, this is the ChIP pipeline. Can I get some summits?')"
	params:
		moduleVer1 = modules['macsVer'],
		moduleVer2 = modules['bedtoolsVer'],
		prefix = '{sample}',
		repeat = REPEATPATH
	shell:
		"""
		module purge && module load {params.moduleVer1} && module load {params.moduleVer2}
		macs2 callpeak -g mm -q 0.01 -f BED -n {params.prefix} -t {input}
		awk 'BEGIN{OFS="\t"}{print $1, $2-50, $3+50, $4, $5, $6}' {output.summits} > {output.exp50}
		bedtools intersect -a {output.exp50} -b {params.repeat} -v -wa > {output.noreps}
		awk '{print $1, $2+50, $3-50, $4, $5, $6}' {output.noreps} > {output.shrink}
		"""

rule makeBedGraph:
	input:
		bed = 'Bed/{sample}-fixall.bed',
		normFactor = 'Counts/{sample}-NormFactor.txt'
	output:
		bedgraph = 'Bedgraph/{sample}-fixall.bedgraph',
		sort = 'Bedgraph/{sample}-sort.bedgraph'
	message: "You're almost to the end of the pipeline! Your bedgraph is being generated, normalized, and sorted..."
	params:
		chrom = chromSizes
	envmodules:
		modules['bedtoolsVer']
	shell:
		"""
		normFactor=$(cat {input.normFactor})
		bedtools genomecov -bga -g {params.chrom} -scale $normFactor -i {input.bed} > {output.bedgraph}
		sort -k1,1 -k2,2n {output.bedgraph} > {output.sort}
		"""

rule convertToBigwig:
	input:
		'Bedgraph/{sample}-sort.bedgraph'
	output:
		'Bigwig/{sample}.bw'
	message: "Congratulations!!! Your sample made it to the end of the pipeline. Your bigwig will be finished shortly..."
	params:
		chrom = chromSizes
	envmodules:
		modules['ucscVer']
	shell:
		"""
		bedGraphToBigWig {input} ${params.chrom} {output}
		"""
