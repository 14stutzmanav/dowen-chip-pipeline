# CUT&RUN Pipeline
## Authors: Megan Justice, Alexis Stutzman

## Quick Start:
#### ☝️ 1. Clone pipeline. (Current stable version: v1.0.0)
```
git clone git@github.com:14stutzmanav/dowen-chip-pipeline.git --branch v1.0.0 --depth 1 && cd dowen-chip-pipeline/ && rm -rf .git
```

Alternatively (if you're in the JDowen lab), you can copy `/proj/dowenlab/projects/new-dowen-ChIP-pipeline` to wherever you're working (probably `pine` for this) and then cd into new-dowen-ChIP-pipeline. Here's how you can do that in a single line:
```
cp -r /proj/dowenlab/projects/new-dowen-ChIP-pipeline ./ && cd new-dowen-ChIP-pipeline
```

#### ✌️ 2. Create `samplesheet.csv` (see below) with descriptive columns of data. 
You can make `samplesheet.csv` in excel and save the file as a `.csv` at the end. Each sample needs two rows of data (one for each end, R1 and R2). The example below (which is tab-deliminated for your ease of reading... note that the `.csv` would turn each tab to a `,`) is set up to run each lane (denoted L001 and L002 in the fastq file from HTSF) separately. 
```
sampleName	fastq	htsfFile
mySample-lane1	htsf_L001_R1.fastq.gz	path/to/htsf_L001_R1.fastq.gz
mySample-lane1	htsf_L001_R2.fastq.gz	path/to/htsf_L001_R2.fastq.gz
mySample-lane2  htsf_L002_R1.fastq.gz	path/to/htsf_L002_R1.fastq.gz
mySample-lane2  htsf_L002_R2.fastq.gz	path/to/htsf_L002_R2.fastq.gz
```
If you want the pipeline to combine multiple fastq files into a single "sample" in the pipeline, just give the fastq files identical sample names. For example, you may want to combine two fastq files that originate from the same ChIP sample but were sequenced across two lanes by HTSF. To combine them in the pipeline, the samplesheet enteries would look something like this:
```
sampleName	fastq	htsfFile
mySample-merged	htsf_L001_R1.fastq.gz	path/to/htsf_L001_R1.fastq.gz
mySample-merged	htsf_L001_R2.fastq.gz	path/to/htsf_L001_R2.fastq.gz
mySample-merged  htsf_L002_R1.fastq.gz	path/to/htsf_L002_R1.fastq.gz
mySample-merged  htsf_L002_R2.fastq.gz	path/to/htsf_L002_R2.fastq.gz
```

#### 🌴 3. Upload and point pipeline to your `samplesheet.csv`.
First, put your `sampleSheet.csv` on longleaf. Once you've done that, you need to edit the `Snakefile.py` in the line where it specifies the path to your sampleSheet. In the stable version of this pipeline, that happens around line 29. Here's the line you're looking to edit (you want to change `example-samplesheet.csv` to whatever your `samplesheet.csv` is named):
```
sampleSheetPath = str('example-samplesheet.csv')
```

#### 🐍 4. Load python.
```
module load python/3.8.8
```

#### 🖐️ 5. Submit the pipeline.
**Note:** If this is your first time using this sampleSheet.csv, it's a good idea to make sure everything is set up properly before you submit your pipeline. You can do so by first running a "dry run."
```
bash dryrunSnakemake.sh Snakefile.py
```
If you see a bunch of text (a lot of it may appear green, depending on your display settings) that ends with a list of jobs and counts, you're good to go! If you didn't something's up. A good place to start debugging is your `sampleSheet.csv` (Questions to ask: are the column names correct? Did I save this as a csv? Are there any extra lines comprised of just commas?).

If you're ready to start the pipeline, you can do so with the following submission:
```
bash runSnakemake.sh Snakefile.py
```

**Alex's Personal Recommendation:** Instead of just submitting the submission script, you should use nohup. This lets you disconnect from longleaf and writes the contents of terminal as jobs are running to an output file called nohup.out. Having nohup.out saves a lot of hassle when you're troubleshooting because it tells you exactly which job caused the error. To submit your pipeline using nohup, run the following command:
```
nohup bash runSnakemake.sh Snakefile.py & disown
```


## Quick Troubleshooting:
#### 🤔 If you're being told `snakmake` isn't found...
...you forgot to load `python`. Check out step 4 then try to run the pipeline again. Good luck! 🤞🍀

#### 🤔 If you're being told there is already an instance of `snakemake` happening in your directory...
1. Make absolutely certain no jobs from this pipeline are running in SLURM.
2. Now, you need to unlock your working directory. To do so, run this:
```
bash unlockSnakemake.sh Snakefile.py
```
Okay, now try to run the pipeline again. Good luck! 🤞🍀

#### 🤔 If you find that a job timed-out...
First, you need to edit `slurmConfig.json` to increase the time allowed for jobs. Here's what that file looks like:
```
{
	"__default__" :
	{
		"time" : "00:30:00",
		"threads" : "1",
		"mem" : "8G"
	},

	"align" :
	{
		"time"	: "10:00:00",
		"threads" : "12",
		"mem" : "100G"
	},

	"convertToBigwig" :
	{
		"mail-type" : "end"
	}

}
```

Each set of bracketed parts with a title corresponds to different rules. At the top, you can see there's a set of paramaters called "__default__". This means that unless the rule is specifically named in a subsequent set of paramaters (so `align` and `convertToBigwig` in this example), then the default paramters will be used. Note that if `align` and `convertToBigwig` don't change the variables listed in `"__default__"`, they will still use the `"__default__"` paramaters.

To increase the time needed for your rule, the best way to do that is setting up a new set of paramaters for that rule. This makes sure that other rules don't run with your increased time, which could penalize you in the SLURM queue and make your pipeline run slower. Note that if you need to inrease the time limit in align, you don't need to add a new section because `slurmConfig.json` already has paramaters listed for that rule. You can just edit the time requirements. Let's say, however, you need to add time to a different rule- as an example, we'll say that we need to let `callPeaks` run for an hour instead of 30min (like it would if it had default parameters). To do that, you would just edit slurmConfig.json to have the following:
```
{
	"__default__" :
	{
		"time" : "00:30:00",
		"threads" : "1",
		"mem" : "8G"
	},

	"align" :
	{
		"time"	: "10:00:00",
		"threads" : "12",
		"mem" : "100G"
	},

	"convertToBigwig" :
	{
		"mail-type" : "end"
	},
  
  	"callPeaks" :
	{
		"time" : "1:00:00"
	}

}
```
Here, I added the paramaters I want for `callPeaks` to the end. Remember to add a comma after the `convertToBigwig` paramaters and make sure the brakets to all slurmConfig paramaters closes after the rule you just added.

Okay, now you're ready to give the pipeline another shot. You'll probably have to unlock your working directory. To do so, make sure there are no jobs from this directory running in SLURM. Then, you unlock by running:
```
bash unlockSnakemake.sh Snakefile.py
```

Okay, now try to submit the pipeline again. Good luck! 🤞🍀
