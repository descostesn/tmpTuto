# An example of what-the-hell error in Snakemake

I want in this post to describe a bug that I got recently and that led me to use constraints on wildcards. The error was not intuitive and can be very disturbing at first if one is not used to the Snakemake logic. Below I will describe step by step how to build a minimal working example that generates the error, what I did to isolate the problem and finally how I solved it.

## Preparing the files

Run the following script to create the folder structure:

```
#!/usr/bin/bash

## Create the project and environment folders
mkdir demo-constraints
cd demo-constraints
mkdir envs

## create the snakeFile that will connect the subworkflows snakeA and snakeB
touch Snakefile
touch config.yaml
```

For more information about the `YAML` format see [here](https://www.cloudbees.com/blog/yaml-tutorial-everything-you-need-get-started). For the next step, you need to have [mamba](https://mamba.readthedocs.io/en/latest/installation.html) installed on your computer. Run the commands:

```
mamba create -p envs/smake
mamba activate envs/smake
mamba install snakemake=7.15.2
```

You can verify that the correct version of snakemake has been installed with:

```
#!/usr/bin/bash

snakemake --version
```


## The workflow

In this workflow, we will code two rules to download fastq files from [Encode](https://www.encodeproject.org/). The URLs will be retrieved from the config file as the information about the species and the associated protocols.

Add the following content to `config.yaml`:

```
testDatasets:
  technique: ["RNASeq", "ChIPSeq", "ATACSeq"]
  organism: ["Mus_musculus", "Homo_sapiens"]
  RNASeq:
    Mus_musculus:
      nameSingleEnd: "limbPolyAPlus"
      singleEnded: "https://www.encodeproject.org/files/ENCFF678XFK/@@download/ENCFF678XFK.fastq.gz"
      namePairedEnd: "CD4PolyAPlus"
      pairedEnded1: "https://www.encodeproject.org/files/ENCFF817RHL/@@download/ENCFF817RHL.fastq.gz"
      pairedEnded2: "https://www.encodeproject.org/files/ENCFF533VQX/@@download/ENCFF533VQX.fastq.gz"
    Homo_sapiens:
      nameSingleEnd: "GM12878PolyAPlus"
      singleEnded: "https://www.encodeproject.org/files/ENCFF729YAX/@@download/ENCFF729YAX.fastq.gz"
      namePairedEnd: "adrenalGlandPolyAPlus"
      pairedEnded1: "https://www.encodeproject.org/files/ENCFF028DUO/@@download/ENCFF028DUO.fastq.gz"
      pairedEnded2: "https://www.encodeproject.org/files/ENCFF470RWW/@@download/ENCFF470RWW.fastq.gz"
  ChIPSeq:
    Mus_musculus:
      nameSingleEnd: "H3K27acMacro"
      singleEnded: "https://www.encodeproject.org/files/ENCFF937IMG/@@download/ENCFF937IMG.fastq.gz"
      namePairedEnd: "H3K27me3Patski"
      pairedEnded1: "https://www.encodeproject.org/files/ENCFF090PQE/@@download/ENCFF090PQE.fastq.gz"
      pairedEnded2: "https://www.encodeproject.org/files/ENCFF362BSH/@@download/ENCFF362BSH.fastq.gz"
    Homo_sapiens:
      nameSingleEnd: "H3K36me3BlaER1"
      singleEnded: "https://www.encodeproject.org/files/ENCFF354RIC/@@download/ENCFF354RIC.fastq.gz"
      namePairedEnd: "BMI1MCF7"
      pairedEnded1: "https://www.encodeproject.org/files/ENCFF825HMN/@@download/ENCFF825HMN.fastq.gz"
      pairedEnded2: "https://www.encodeproject.org/files/ENCFF240ZBS/@@download/ENCFF240ZBS.fastq.gz"
  ATACSeq:
    Mus_musculus:
      nameSingleEnd: "ATACErythroblast"
      singleEnded: "https://www.encodeproject.org/files/ENCFF535AXY/@@download/ENCFF535AXY.fastq.gz"
      namePairedEnd: "ATACPatski"
      pairedEnded1: "https://www.encodeproject.org/files/ENCFF600AUS/@@download/ENCFF600AUS.fastq.gz"
      pairedEnded2: "https://www.encodeproject.org/files/ENCFF273QXE/@@download/ENCFF273QXE.fastq.gz"
    Homo_sapiens:
      nameSingleEnd: "ATACA549"
      singleEnded: "https://www.encodeproject.org/files/ENCFF022FVS/@@download/ENCFF022FVS.fastq.gz"
      namePairedEnd: "ATACCD4positive"
      pairedEnded1: "https://www.encodeproject.org/files/ENCFF526FOQ/@@download/ENCFF526FOQ.fastq.gz"
      pairedEnded2: "https://www.encodeproject.org/files/ENCFF310QJA/@@download/ENCFF310QJA.fastq.gz"
```

Define the workflow in `Snakefile`:

```
import pandas as pd

configfile: "config.yaml"

onstart:
    print("##### DOWNLOAD FASTQ FIILES #####\n") 


###############################################################################
# Creating input table
###############################################################################

# Build the table of test datasets to download
samplesData = []

for tech in config["testDatasets"]["technique"]:
  for org in config["testDatasets"]["organism"]:
    pathSingle = config["testDatasets"][tech][org]["singleEnded"]
    nameSingle = config["testDatasets"][tech][org]["nameSingleEnd"]
    pathPaired1 = config["testDatasets"][tech][org]["pairedEnded1"]
    pathPaired2 = config["testDatasets"][tech][org]["pairedEnded2"]
    namePaired = config["testDatasets"][tech][org]["namePairedEnd"]
    samplesData.append([nameSingle, tech, org, "single", pathSingle, "NA"])
    samplesData.append([namePaired, tech, org, "paired", pathPaired1, pathPaired2])

df = pd.DataFrame(samplesData)
df.rename(columns={0: 'samples', 1: 'library_strategy', 2: 'organism', 3: 'library_layout', 4: 'link1', 5: 'link2'}, inplace=True)


###############################################################################
# Variables definition
###############################################################################

# Splitting the table into single or paired end experiments

index_single = df['library_layout'] == 'single'
index_paired = df['library_layout'] == 'paired'
df_single = df[index_single]
df_paired = df[index_paired]

# Output files names

SINGLESAMPLES = df_single['samples'].tolist()
PAIREDSAMPLES = df_paired['samples'].tolist()

# For Retrieving links to download sra files

samples_single_forlinks = pd.DataFrame(df_single).set_index("samples",drop=False)
samples_paired_forlinks = pd.DataFrame(df_paired).set_index("samples",drop=False)

# Technique names
SINGLETECH = df_single['library_strategy'].tolist()
PAIREDTECH = df_paired['library_strategy'].tolist()

## Species name
SPECIESSINGLE = df_single['organism'].tolist()
SPECIESPAIRED = df_paired['organism'].tolist()

## Layout names
LAYOUTSINGLE = df_single['library_layout'].tolist()
LAYOUTPAIRED = df_paired['library_layout'].tolist()


############
# Rule all
############

rule all:
  input:
    expand("results/{speciessingle}/fastq/{techniquesingle}/{layoutsingle}/allchrom/{samplenamesingle}.fastq.gz", zip, speciessingle=SPECIESSINGLE, techniquesingle=SINGLETECH, layoutsingle=LAYOUTSINGLE, samplenamesingle=SINGLESAMPLES),
    expand("results/{speciespaired}/fastq/{techniquepaired}/{layoutpaired}/allchrom/{samplenamepaired}_1.fastq.gz", zip, speciespaired=SPECIESPAIRED, techniquepaired=PAIREDTECH, layoutpaired=LAYOUTPAIRED, samplenamepaired=PAIREDSAMPLES),
    expand("results/{speciespaired}/fastq/{techniquepaired}/{layoutpaired}/allchrom/{samplenamepaired}_2.fastq.gz", zip, speciespaired=SPECIESPAIRED, techniquepaired=PAIREDTECH, layoutpaired=LAYOUTPAIRED, samplenamepaired=PAIREDSAMPLES)


rule download_fastq_single:
  output:
    "results/{speciessingle}/fastq/{techniquesingle}/{layoutsingle}/allchrom/{samplenamesingle}.fastq.gz"
  params:
    outputdirectory = lambda wildcards: f"results/{wildcards.speciessingle}/fastq/{wildcards.techniquesingle}/{wildcards.layoutsingle}/fastq/allchrom",
    linksingle = lambda wildcards: samples_single_forlinks.loc[wildcards.samplenamesingle, "link1"]
  threads: 1
  shell:
    """
    echo "Downloading {params.linksingle}"
    wget --directory-prefix={params.outputdirectory} {params.linksingle}
    sleep 10s
    FILENAME=`basename {params.linksingle}`
    mv {params.outputdirectory}/$FILENAME {output}  
    """

rule download_fastq_paired:
  output:
    pair1 = "results/{speciespaired}/fastq/{techniquepaired}/{layoutpaired}/allchrom/{samplenamepaired}_1.fastq.gz",
    pair2 = "results/{speciespaired}/fastq/{techniquepaired}/{layoutpaired}/allchrom/{samplenamepaired}_2.fastq.gz"
  threads: 1
  params:
    outputdirectory = lambda wildcards: f"results/{wildcards.speciespaired}/fastq/{wildcards.techniquepaired}/{wildcards.layoutpaired}/fastq/allchrom",
    linkpair1 = lambda wildcards: samples_paired_forlinks.loc[wildcards.samplenamepaired, "link1"],
    linkpair2 = lambda wildcards: samples_paired_forlinks.loc[wildcards.samplenamepaired, "link2"]
  shell:
    """
    echo "Downloading {params.linkpair1} and {params.linkpair2}"
    wget --directory-prefix={params.outputdirectory} {params.linkpair1}
    wget --directory-prefix={params.outputdirectory} {params.linkpair2}
    sleep 10s
    FILENAME1=`basename {params.linkpair1}`
    FILENAME2=`basename {params.linkpair2}`
    mv {params.outputdirectory}/$FILENAME1 {output.pair1}
    mv {params.outputdirectory}/$FILENAME2 {output.pair2}
    """
```

## The workflow explained step by step

The Snakefile starts with:

```
import pandas as pd

configfile: "config.yaml"

onstart:
    print("##### DOWNLOAD FASTQ FIILES #####\n") 
```

These commands first import the python module `pandas` that will be used with the prefix `pd`. The path to the configuration file is then indicated with `configfile: "config.yaml"`. The configuration file (see above) contains information about the different techniques that were used to generate the data (`technique: ["RNASeq", "ChIPSeq", "ATACSeq"]`), the organism from which the data were generated (`organism: ["Mus_musculus", "Homo_sapiens"]`) and the URLs are further organised according to the techniques and organisms. The `onstart` section prints a message in the terminal upon workflow invocation.


The `Creating input table` section builds a panda dataframe from the `config.yaml` with the following information:

```
                  samples library_strategy      organism library_layout                                              link1                                              link2
0           limbPolyAPlus           RNASeq  Mus_musculus         single  https://www.encodeproject.org/files/ENCFF678XF...                                                 NA
1            CD4PolyAPlus           RNASeq  Mus_musculus         paired  https://www.encodeproject.org/files/ENCFF817RH...  https://www.encodeproject.org/files/ENCFF533VQ...
2        GM12878PolyAPlus           RNASeq  Homo_sapiens         single  https://www.encodeproject.org/files/ENCFF729YA...                                                 NA
3   adrenalGlandPolyAPlus           RNASeq  Homo_sapiens         paired  https://www.encodeproject.org/files/ENCFF028DU...  https://www.encodeproject.org/files/ENCFF470RW...
4            H3K27acMacro          ChIPSeq  Mus_musculus         single  https://www.encodeproject.org/files/ENCFF937IM...                                                 NA
5          H3K27me3Patski          ChIPSeq  Mus_musculus         paired  https://www.encodeproject.org/files/ENCFF090PQ...  https://www.encodeproject.org/files/ENCFF362BS...
6          H3K36me3BlaER1          ChIPSeq  Homo_sapiens         single  https://www.encodeproject.org/files/ENCFF354RI...                                                 NA
7                BMI1MCF7          ChIPSeq  Homo_sapiens         paired  https://www.encodeproject.org/files/ENCFF825HM...  https://www.encodeproject.org/files/ENCFF240ZB...
8        ATACErythroblast          ATACSeq  Mus_musculus         single  https://www.encodeproject.org/files/ENCFF535AX...                                                 NA
9              ATACPatski          ATACSeq  Mus_musculus         paired  https://www.encodeproject.org/files/ENCFF600AU...  https://www.encodeproject.org/files/ENCFF273QX...
10               ATACA549          ATACSeq  Homo_sapiens         single  https://www.encodeproject.org/files/ENCFF022FV...                                                 NA
11        ATACCD4positive          ATACSeq  Homo_sapiens         paired  https://www.encodeproject.org/files/ENCFF526FO...  https://www.encodeproject.org/files/ENCFF310QJ...
```

The rows are then separated into two tables of `single` and `paired` samples with:

```
# Splitting the table into single or paired end experiments

index_single = df['library_layout'] == 'single'
index_paired = df['library_layout'] == 'paired'
df_single = df[index_single]
df_paired = df[index_paired]
```

The data frames are then indexed in order to be able to retrieve the URLs with the samples names (`sample`) and the different variables used to build the output paths are created:

```
# Output files names

SINGLESAMPLES = df_single['samples'].tolist()
PAIREDSAMPLES = df_paired['samples'].tolist()

# For Retrieving links to download sra files

samples_single_forlinks = pd.DataFrame(df_single).set_index("samples",drop=False)
samples_paired_forlinks = pd.DataFrame(df_paired).set_index("samples",drop=False)

# Technique names
SINGLETECH = df_single['library_strategy'].tolist()
PAIREDTECH = df_paired['library_strategy'].tolist()

## Species name
SPECIESSINGLE = df_single['organism'].tolist()
SPECIESPAIRED = df_paired['organism'].tolist()

## Layout names
LAYOUTSINGLE = df_single['library_layout'].tolist()
LAYOUTPAIRED = df_paired['library_layout'].tolist()
```

Finally, `rule all` contains the files to be created and the rules `download_fastq_single`/`download_fastq_paired` use the indexed data frames `samples_single_forlinks`/`samples_paired_forlinks` to download the fastq files.



## Running the pipeline and error


To dry-run the pipeline use:

```
#!/usr/bin/bash

snakemake --cores 1 -n
```


You should see the error:

```
AmbiguousRuleException:
Rules download_fastq_paired and download_fastq_single are ambiguous for the file results/Mus_musculus/fastq/RNASeq/paired/allchrom/CD4PolyAPlus_1.fastq.gz.
Consider starting rule output with a unique prefix, constrain your wildcards, or use the ruleorder directive.
Wildcards:
    download_fastq_paired: layoutpaired=paired,samplenamepaired=CD4PolyAPlus,speciespaired=Mus_musculus,techniquepaired=RNASeq
    download_fastq_single: layoutsingle=paired,samplenamesingle=CD4PolyAPlus_1,speciessingle=Mus_musculus,techniquesingle=RNASeq
Expected input files:
    download_fastq_paired: 
    download_fastq_single: Expected output files:
    download_fastq_paired: results/Mus_musculus/fastq/RNASeq/paired/allchrom/CD4PolyAPlus_1.fastq.gz results/Mus_musculus/fastq/RNASeq/paired/allchrom/CD4PolyAPlus_2.fastq.gz
    download_fastq_single: results/Mus_musculus/fastq/RNASeq/paired/allchrom/CD4PolyAPlus_1.fastq.gz
```

The first thing to notice is that both rules, even if built on two different dataframes, process the same sample (`results/Mus_musculus/fastq/RNASeq/paired/allchrom/CD4PolyAPlus_1.fastq.gz`). Looking at the `Wildcards` section, you will notice that the rule dedicated to the single-ended samples is now receiving `paired` parameters (`layoutsingle=paired`). More interestingly, the `download_fastq_single` rule does not contain the correct wildcards `samplenamesingle` and the `_1` suffix was added. This suffix is only present in the output of `download_fastq_paired`. Notice also that the wildcards `samplenamepaired=CD4PolyAPlus` is correct as the expected output files (`results/Mus_musculus/fastq/RNASeq/paired/allchrom/CD4PolyAPlus_1.fastq.gz`/`results/Mus_musculus/fastq/RNASeq/paired/allchrom/CD4PolyAPlus_2.fastq.gz`).


=> The question is how can the paired samples could be injected to the `download_fastq_single`? Why does the '_1' suffix suddenly appeared in that rule?


## Simplifying the workflow to isolate the problem


Remove everything related to the paired rule to generate only the files of `download_fastq_single`. Then simplify the rule to a minimal working example. Copy the following content into a `Snakefile2`:

```
import pandas as pd

configfile: "config.yaml"

onstart:
    print("##### DOWNLOAD FASTQ FIILES #####\n") 


###############################################################################
# Creating input table
###############################################################################

# Build the table of test datasets to download
samplesData = []

for tech in config["testDatasets"]["technique"]:
  for org in config["testDatasets"]["organism"]:
    pathSingle = config["testDatasets"][tech][org]["singleEnded"]
    nameSingle = config["testDatasets"][tech][org]["nameSingleEnd"]
    pathPaired1 = config["testDatasets"][tech][org]["pairedEnded1"]
    pathPaired2 = config["testDatasets"][tech][org]["pairedEnded2"]
    namePaired = config["testDatasets"][tech][org]["namePairedEnd"]
    samplesData.append([nameSingle, tech, org, "single", pathSingle, "NA"])
    samplesData.append([namePaired, tech, org, "paired", pathPaired1, pathPaired2])

df = pd.DataFrame(samplesData)
df.rename(columns={0: 'samples', 1: 'library_strategy', 2: 'organism', 3: 'library_layout', 4: 'link1', 5: 'link2'}, inplace=True)


###############################################################################
# Variables definition
###############################################################################

# Splitting the table into single or paired end experiments
index_single = df['library_layout'] == 'single'
df_single = df[index_single]

# Output files names
SINGLESAMPLES = df_single['samples'].tolist()

# For Retrieving links to download sra files
samples_single_forlinks = pd.DataFrame(df_single).set_index("samples",drop=False)

# Technique names
SINGLETECH = df_single['library_strategy'].tolist()

## Species name
SPECIESSINGLE = df_single['organism'].tolist()

## Layout names
LAYOUTSINGLE = df_single['library_layout'].tolist()

############
# Rule all
############

rule all:
  input:
    expand("results/{speciessingle}/fastq/{techniquesingle}/{layoutsingle}/allchrom/{samplenamesingle}.fastq.gz", zip, speciessingle=SPECIESSINGLE, techniquesingle=SINGLETECH, layoutsingle=LAYOUTSINGLE, samplenamesingle=SINGLESAMPLES)

rule download_fastq_single:
  output:
    "results/{speciessingle}/fastq/{techniquesingle}/{layoutsingle}/allchrom/{samplenamesingle}.fastq.gz"
  threads: 1
  shell:
    "echo 'hello' > {output}"
```

Perform a dry-run:

```
#!/usr/bin/bash

snakemake --cores 1 --snakefile Snakefile2 -n
```

You should obtain:

```
Building DAG of jobs...
Job stats:
job                      count    min threads    max threads
---------------------  -------  -------------  -------------
all                          1              1              1
download_fastq_single        6              1              1
total                        7              1              1


[Thu Oct 13 15:46:16 2022]
rule download_fastq_single:
    output: results/Homo_sapiens/fastq/RNASeq/single/allchrom/GM12878PolyAPlus.fastq.gz
    jobid: 2
    reason: Missing output files: results/Homo_sapiens/fastq/RNASeq/single/allchrom/GM12878PolyAPlus.fastq.gz
    wildcards: speciessingle=Homo_sapiens, techniquesingle=RNASeq, layoutsingle=single, samplenamesingle=GM12878PolyAPlus
    resources: tmpdir=/tmp


[Thu Oct 13 15:46:16 2022]
rule download_fastq_single:
    output: results/Mus_musculus/fastq/ChIPSeq/single/allchrom/H3K27acMacro.fastq.gz
    jobid: 3
    reason: Missing output files: results/Mus_musculus/fastq/ChIPSeq/single/allchrom/H3K27acMacro.fastq.gz
    wildcards: speciessingle=Mus_musculus, techniquesingle=ChIPSeq, layoutsingle=single, samplenamesingle=H3K27acMacro
    resources: tmpdir=/tmp


[Thu Oct 13 15:46:16 2022]
rule download_fastq_single:
    output: results/Homo_sapiens/fastq/ChIPSeq/single/allchrom/H3K36me3BlaER1.fastq.gz
    jobid: 4
    reason: Missing output files: results/Homo_sapiens/fastq/ChIPSeq/single/allchrom/H3K36me3BlaER1.fastq.gz
    wildcards: speciessingle=Homo_sapiens, techniquesingle=ChIPSeq, layoutsingle=single, samplenamesingle=H3K36me3BlaER1
    resources: tmpdir=/tmp


[Thu Oct 13 15:46:16 2022]
rule download_fastq_single:
    output: results/Mus_musculus/fastq/ATACSeq/single/allchrom/ATACErythroblast.fastq.gz
    jobid: 5
    reason: Missing output files: results/Mus_musculus/fastq/ATACSeq/single/allchrom/ATACErythroblast.fastq.gz
    wildcards: speciessingle=Mus_musculus, techniquesingle=ATACSeq, layoutsingle=single, samplenamesingle=ATACErythroblast
    resources: tmpdir=/tmp


[Thu Oct 13 15:46:16 2022]
rule download_fastq_single:
    output: results/Homo_sapiens/fastq/ATACSeq/single/allchrom/ATACA549.fastq.gz
    jobid: 6
    reason: Missing output files: results/Homo_sapiens/fastq/ATACSeq/single/allchrom/ATACA549.fastq.gz
    wildcards: speciessingle=Homo_sapiens, techniquesingle=ATACSeq, layoutsingle=single, samplenamesingle=ATACA549
    resources: tmpdir=/tmp


[Thu Oct 13 15:46:16 2022]
rule download_fastq_single:
    output: results/Mus_musculus/fastq/RNASeq/single/allchrom/limbPolyAPlus.fastq.gz
    jobid: 1
    reason: Missing output files: results/Mus_musculus/fastq/RNASeq/single/allchrom/limbPolyAPlus.fastq.gz
    wildcards: speciessingle=Mus_musculus, techniquesingle=RNASeq, layoutsingle=single, samplenamesingle=limbPolyAPlus
    resources: tmpdir=/tmp


[Thu Oct 13 15:46:16 2022]
localrule all:
    input: results/Mus_musculus/fastq/RNASeq/single/allchrom/limbPolyAPlus.fastq.gz, results/Homo_sapiens/fastq/RNASeq/single/allchrom/GM12878PolyAPlus.fastq.gz, results/Mus_musculus/fastq/ChIPSeq/single/allchrom/H3K27acMacro.fastq.gz, results/Homo_sapiens/fastq/ChIPSeq/single/allchrom/H3K36me3BlaER1.fastq.gz, results/Mus_musculus/fastq/ATACSeq/single/allchrom/ATACErythroblast.fastq.gz, results/Homo_sapiens/fastq/ATACSeq/single/allchrom/ATACA549.fastq.gz
    jobid: 0
    reason: Input files updated by another job: results/Mus_musculus/fastq/RNASeq/single/allchrom/limbPolyAPlus.fastq.gz, results/Homo_sapiens/fastq/ATACSeq/single/allchrom/ATACA549.fastq.gz, results/Homo_sapiens/fastq/RNASeq/single/allchrom/GM12878PolyAPlus.fastq.gz, results/Homo_sapiens/fastq/ChIPSeq/single/allchrom/H3K36me3BlaER1.fastq.gz, results/Mus_musculus/fastq/ChIPSeq/single/allchrom/H3K27acMacro.fastq.gz, results/Mus_musculus/fastq/ATACSeq/single/allchrom/ATACErythroblast.fastq.gz
    resources: tmpdir=/tmp

Job stats:
job                      count    min threads    max threads
---------------------  -------  -------------  -------------
all                          1              1              1
download_fastq_single        6              1              1
total                        7              1              1

Reasons:
    (check individual jobs above for details)
    input files updated by another job:
        all
    missing output files:
        download_fastq_single

This was a dry-run (flag -n). The order of jobs does not reflect the order of execution.
```

You can see above that the `single` samples are processed correctly. We can deduce that the problem was coming from the `download_fastq_paired`. Let's now add back the paired rule but with a minimal structure to reduce the number of possible sources of error. We want first to see if the error could come from the `output` section of the rule. Copy the following content into `Snakefile3`:

```
import pandas as pd

configfile: "config.yaml"

onstart:
    print("##### DOWNLOAD FASTQ FIILES #####\n") 


###############################################################################
# Creating input table
###############################################################################

# Build the table of test datasets to download
samplesData = []

for tech in config["testDatasets"]["technique"]:
  for org in config["testDatasets"]["organism"]:
    pathSingle = config["testDatasets"][tech][org]["singleEnded"]
    nameSingle = config["testDatasets"][tech][org]["nameSingleEnd"]
    pathPaired1 = config["testDatasets"][tech][org]["pairedEnded1"]
    pathPaired2 = config["testDatasets"][tech][org]["pairedEnded2"]
    namePaired = config["testDatasets"][tech][org]["namePairedEnd"]
    samplesData.append([nameSingle, tech, org, "single", pathSingle, "NA"])
    samplesData.append([namePaired, tech, org, "paired", pathPaired1, pathPaired2])

df = pd.DataFrame(samplesData)
df.rename(columns={0: 'samples', 1: 'library_strategy', 2: 'organism', 3: 'library_layout', 4: 'link1', 5: 'link2'}, inplace=True)


###############################################################################
# Variables definition
###############################################################################

# Splitting the table into single or paired end experiments

index_single = df['library_layout'] == 'single'
index_paired = df['library_layout'] == 'paired'
df_single = df[index_single]
df_paired = df[index_paired]

# Output files names

SINGLESAMPLES = df_single['samples'].tolist()
PAIREDSAMPLES = df_paired['samples'].tolist()

# For Retrieving links to download sra files

samples_single_forlinks = pd.DataFrame(df_single).set_index("samples",drop=False)
samples_paired_forlinks = pd.DataFrame(df_paired).set_index("samples",drop=False)

# Technique names
SINGLETECH = df_single['library_strategy'].tolist()
PAIREDTECH = df_paired['library_strategy'].tolist()

## Species name
SPECIESSINGLE = df_single['organism'].tolist()
SPECIESPAIRED = df_paired['organism'].tolist()

## Layout names
LAYOUTSINGLE = df_single['library_layout'].tolist()
LAYOUTPAIRED = df_paired['library_layout'].tolist()


############
# Rule all
############

rule all:
  input:
    expand("results/{speciessingle}/fastq/{techniquesingle}/{layoutsingle}/allchrom/{samplenamesingle}.fastq.gz", zip, speciessingle=SPECIESSINGLE, techniquesingle=SINGLETECH, layoutsingle=LAYOUTSINGLE, samplenamesingle=SINGLESAMPLES),
    expand("results/{speciespaired}/fastq/{techniquepaired}/{layoutpaired}/allchrom/{samplenamepaired}_1.fastq.gz", zip, speciespaired=SPECIESPAIRED, techniquepaired=PAIREDTECH, layoutpaired=LAYOUTPAIRED, samplenamepaired=PAIREDSAMPLES)


rule download_fastq_single:
  output:
    "results/{speciessingle}/fastq/{techniquesingle}/{layoutsingle}/allchrom/{samplenamesingle}.fastq.gz"
  threads: 1
  shell:
    "echo 'hello' > {output}"

rule download_fastq_paired:
  output:
    "results/{speciespaired}/fastq/{techniquepaired}/{layoutpaired}/allchrom/{samplenamepaired}_1.fastq.gz"
  threads: 1
  shell:
    "echo 'hello' > {output}"
```

Perform a dry-run:

```
#!/usr/bin/bash

snakemake --cores 1 --snakefile Snakefile3 -n
```

You can see that the error is back:

```
Building DAG of jobs...
AmbiguousRuleException:
Rules download_fastq_paired and download_fastq_single are ambiguous for the file results/Mus_musculus/fastq/RNASeq/paired/allchrom/CD4PolyAPlus_1.fastq.gz.
Consider starting rule output with a unique prefix, constrain your wildcards, or use the ruleorder directive.
Wildcards:
    download_fastq_paired: layoutpaired=paired,samplenamepaired=CD4PolyAPlus,speciespaired=Mus_musculus,techniquepaired=RNASeq
    download_fastq_single: layoutsingle=paired,samplenamesingle=CD4PolyAPlus_1,speciessingle=Mus_musculus,techniquesingle=RNASeq
Expected input files:
    download_fastq_paired: 
    download_fastq_single: Expected output files:
    download_fastq_paired: results/Mus_musculus/fastq/RNASeq/paired/allchrom/CD4PolyAPlus_1.fastq.gz
    download_fastq_single: results/Mus_musculus/fastq/RNASeq/paired/allchrom/CD4PolyAPlus_1.fastq.gz
```

We can now be confident that the problem comes from the line `"results/{speciespaired}/fastq/{techniquepaired}/{layoutpaired}/allchrom/{samplenamepaired}_1.fastq.gz"`. To be even more sure about it, let's change this line (use a `Snakefile4`):

```
import pandas as pd

configfile: "config.yaml"

onstart:
    print("##### DOWNLOAD FASTQ FIILES #####\n") 


###############################################################################
# Creating input table
###############################################################################

# Build the table of test datasets to download
samplesData = []

for tech in config["testDatasets"]["technique"]:
  for org in config["testDatasets"]["organism"]:
    pathSingle = config["testDatasets"][tech][org]["singleEnded"]
    nameSingle = config["testDatasets"][tech][org]["nameSingleEnd"]
    pathPaired1 = config["testDatasets"][tech][org]["pairedEnded1"]
    pathPaired2 = config["testDatasets"][tech][org]["pairedEnded2"]
    namePaired = config["testDatasets"][tech][org]["namePairedEnd"]
    samplesData.append([nameSingle, tech, org, "single", pathSingle, "NA"])
    samplesData.append([namePaired, tech, org, "paired", pathPaired1, pathPaired2])

df = pd.DataFrame(samplesData)
df.rename(columns={0: 'samples', 1: 'library_strategy', 2: 'organism', 3: 'library_layout', 4: 'link1', 5: 'link2'}, inplace=True)


###############################################################################
# Variables definition
###############################################################################

# Splitting the table into single or paired end experiments

index_single = df['library_layout'] == 'single'
index_paired = df['library_layout'] == 'paired'
df_single = df[index_single]
df_paired = df[index_paired]

# Output files names

SINGLESAMPLES = df_single['samples'].tolist()
PAIREDSAMPLES = df_paired['samples'].tolist()

# For Retrieving links to download sra files

samples_single_forlinks = pd.DataFrame(df_single).set_index("samples",drop=False)
samples_paired_forlinks = pd.DataFrame(df_paired).set_index("samples",drop=False)

# Technique names
SINGLETECH = df_single['library_strategy'].tolist()
PAIREDTECH = df_paired['library_strategy'].tolist()

## Species name
SPECIESSINGLE = df_single['organism'].tolist()
SPECIESPAIRED = df_paired['organism'].tolist()

## Layout names
LAYOUTSINGLE = df_single['library_layout'].tolist()
LAYOUTPAIRED = df_paired['library_layout'].tolist()


############
# Rule all
############

rule all:
  input:
    expand("results/{speciessingle}/fastq/{techniquesingle}/{layoutsingle}/allchrom/{samplenamesingle}.fastq.gz", zip, speciessingle=SPECIESSINGLE, techniquesingle=SINGLETECH, layoutsingle=LAYOUTSINGLE, samplenamesingle=SINGLESAMPLES),
    "results/test/test.txt"


rule download_fastq_single:
  output:
    "results/{speciessingle}/fastq/{techniquesingle}/{layoutsingle}/allchrom/{samplenamesingle}.fastq.gz"
  threads: 1
  shell:
    "echo 'hello' > {output}"

rule download_fastq_paired:
  output:
    "results/test/test.txt"
  threads: 1
  shell:
    "echo 'hello' > {output}"
```

Perform a dry-run (`snakemake --cores 1 --snakefile Snakefile4 -n`) and you will notice that the error is gone. We can now be sure that the problem is coming from the output section of `download_fastq_paired` and more precisely from its wildcards. In order to isolate which wildcards is problematic, let's test separately the ones from the path and the one from the file name:

Copy this content into `Snakefile5` and perform a dry-run (`snakemake --cores 1 --snakefile Snakefile5 -n`):

```
import pandas as pd

configfile: "config.yaml"

onstart:
    print("##### DOWNLOAD FASTQ FIILES #####\n") 


###############################################################################
# Creating input table
###############################################################################

# Build the table of test datasets to download
samplesData = []

for tech in config["testDatasets"]["technique"]:
  for org in config["testDatasets"]["organism"]:
    pathSingle = config["testDatasets"][tech][org]["singleEnded"]
    nameSingle = config["testDatasets"][tech][org]["nameSingleEnd"]
    pathPaired1 = config["testDatasets"][tech][org]["pairedEnded1"]
    pathPaired2 = config["testDatasets"][tech][org]["pairedEnded2"]
    namePaired = config["testDatasets"][tech][org]["namePairedEnd"]
    samplesData.append([nameSingle, tech, org, "single", pathSingle, "NA"])
    samplesData.append([namePaired, tech, org, "paired", pathPaired1, pathPaired2])

df = pd.DataFrame(samplesData)
df.rename(columns={0: 'samples', 1: 'library_strategy', 2: 'organism', 3: 'library_layout', 4: 'link1', 5: 'link2'}, inplace=True)


###############################################################################
# Variables definition
###############################################################################

# Splitting the table into single or paired end experiments

index_single = df['library_layout'] == 'single'
index_paired = df['library_layout'] == 'paired'
df_single = df[index_single]
df_paired = df[index_paired]

# Output files names

SINGLESAMPLES = df_single['samples'].tolist()
PAIREDSAMPLES = df_paired['samples'].tolist()

# For Retrieving links to download sra files

samples_single_forlinks = pd.DataFrame(df_single).set_index("samples",drop=False)
samples_paired_forlinks = pd.DataFrame(df_paired).set_index("samples",drop=False)

# Technique names
SINGLETECH = df_single['library_strategy'].tolist()
PAIREDTECH = df_paired['library_strategy'].tolist()

## Species name
SPECIESSINGLE = df_single['organism'].tolist()
SPECIESPAIRED = df_paired['organism'].tolist()

## Layout names
LAYOUTSINGLE = df_single['library_layout'].tolist()
LAYOUTPAIRED = df_paired['library_layout'].tolist()


############
# Rule all
############

rule all:
  input:
    expand("results/{speciessingle}/fastq/{techniquesingle}/{layoutsingle}/allchrom/{samplenamesingle}.fastq.gz", zip, speciessingle=SPECIESSINGLE, techniquesingle=SINGLETECH, layoutsingle=LAYOUTSINGLE, samplenamesingle=SINGLESAMPLES),
    expand("results/{speciespaired}/fastq/{techniquepaired}/{layoutpaired}/allchrom/test_1.txt", zip, speciespaired=SPECIESPAIRED, techniquepaired=PAIREDTECH, layoutpaired=LAYOUTPAIRED, samplenamepaired=PAIREDSAMPLES)


rule download_fastq_single:
  output:
    "results/{speciessingle}/fastq/{techniquesingle}/{layoutsingle}/allchrom/{samplenamesingle}.fastq.gz"
  threads: 1
  shell:
    "echo 'hello' > {output}"

rule download_fastq_paired:
  output:
    "results/{speciespaired}/fastq/{techniquepaired}/{layoutpaired}/allchrom/test_1.txt"
  threads: 1
  shell:
    "echo 'hello' > {output}"
```

Since the error was not generated, we can think that the problem is coming from the file name. Let's change `test_1.txt` to `test_1.fastq.gz` (`cp Snakefile5 Snakefile6`, modify the file and run `snakemake --cores 1 --snakefile Snakefile6 -n`). You should obtain:

```
AmbiguousRuleException:
Rules download_fastq_paired and download_fastq_single are ambiguous for the file results/Mus_musculus/fastq/RNASeq/paired/allchrom/test_1.fastq.gz.
Consider starting rule output with a unique prefix, constrain your wildcards, or use the ruleorder directive.
Wildcards:
    download_fastq_paired: layoutpaired=paired,speciespaired=Mus_musculus,techniquepaired=RNASeq
    download_fastq_single: layoutsingle=paired,samplenamesingle=test_1,speciessingle=Mus_musculus,techniquesingle=RNASeq
Expected input files:
    download_fastq_paired: 
    download_fastq_single: Expected output files:
    download_fastq_paired: results/Mus_musculus/fastq/RNASeq/paired/allchrom/test_1.fastq.gz
    download_fastq_single: results/Mus_musculus/fastq/RNASeq/paired/allchrom/test_1.fastq.gz
```

## Interpretation

You can see that the problem observed before is back. That means that when two files have the same structure (`test_1.fastq.gz` and `{samplenamesingle}.fastq.gz`), the back-propagation system of Snakemake can create ambiguity about the files to produce. Where one could see two different structures in:

1) "results/{speciessingle}/fastq/{techniquesingle}/{layoutsingle}/allchrom/{samplenamesingle}.fastq.gz"
2) "results/{speciespaired}/fastq/{techniquepaired}/{layoutpaired}/allchrom/{samplenamepaired}_1.fastq.gz"

Snakemake seems to interpret them as the single pattern: `folder1/folder2/folder3/folder4/folder5/folder6/filename.fastq.gz`. It will then look if this pattern found in `download_fastq_paired` can be found in another rule (here `download_fastq_single`) and will back-propagate the wildcards `filename` to this other rule.


## Fixing the problem

You might have noticed in the `AmbiguousRuleException` the message `Consider starting rule output with a unique prefix, constrain your wildcards, or use the ruleorder directive`. Snakemake offers a way to [constrain wildcards](download_fastq_paired) with the keyword `wildcard_constraints`. The trick here is to remove the possibility of the underscore being included in the wildcards with a regular expression. Add the following code to each rule of `Snakefile` and perform a dry-run (be careful about the indentation when you paste the code below):

```
[...]

rule download_fastq_single:
[...]
  wildcard_constraints:
    samplenamesingle="[0-9A-Za-z]+"
[...]

rule download_fastq_paired:
[...]
  wildcard_constraints:
    samplenamepaired="[0-9A-Za-z]+"
[...]
```

## Conclusion

Having in mind the fact that Snakemake will look at the targets first and then go backward to determine the rules that created them is key to understand some of the errors that one could encounter. We have seen here that the origin of the problem impacting the `download_fastq_single` was actually a downstream rule. Our example being limited to two rules, it was not very complicated to figure it out. This can become really difficult when working with a huge number of rules interconnected in multiple ways. I would say that the key in such a case is to sequentially simplify the DAG to be able to spot the problematic rule.
