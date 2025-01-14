# PGAP2

Pan-Genome Analysis Pipeline 2


## Quick start
### <font style="color:rgb(51, 51, 51);">Basic usage</font>
The input directory contains all the genome and annotation files.

PGAP2 supports multiple input formats: GFF files in the same format as those output by Prokka, GFF files with their corresponding genome FASTA files in separate files, GenBank flat files (GBFF), or just genome FASTA files (with `--reannot` required).

Different formats of input files can be mixed in one input directory. PGAP2 will recognize and process them based on their prefixes and suffixes.

```bash
pgap2 main -i inputdir/ -o outputdir/
```

### Preprocessing
Quality checks and visualization are conducted by PGAP2 during the preprocessing step. PGAP2 generates an interactive HTML file and corresponding vector figures to help users understand their input data. The input data and pre-alignment results are stored as a pickle file for quick restarting of the same calculation step.

```bash
pgap2 prep -i inputdir/ -o outputdir/
```

### Postprocessing
The postprocessing pipeline is performed by PGAP2. There are various submodules integrated into the postprocessing module, such as statistical analysis, single-copy tree building, population clustering, and Tajima's D test. Regardless of which submodule you want to use, you can always run it as follows:

```bash
pgap2 post [submodule] [options] -i inputdir/ -o outputdir/
```

The inputdir is the outputdir of main module.

PGAP2 also support statistical analysis using a PAV file indepandently:

```bash
pgap2 post profile --pav your_pav_file -o outputdir/
```

## Installation
The best way to install full version of PGAP2 package is using conda:

```bash
conda create -n pgap2 -c conda-forge -c bioconda -c defaults pgap2
```

<font style="color:rgb(52, 73, 94);">alternatively it is often faster to use the </font>[mamba](https://github.com/mamba-org/mamba)<font style="color:rgb(52, 73, 94);"> solver:</font>

```bash
conda create -n pgap2 -c conda-forge mamba
conda activate pgap2 
mamba install -c conda-forge -c bioconda -c defaults pgap2
```

Or sometimes you only want to carry out a specific function, such as partioning and don't want install too many extra softwares for fully version of PGAP2, then you can just install PGAP2:

```bash
pip install pgap2
```

Or via source file:

```bash
git clone https://github.com/bucongfan/PGAP2
```

And then install extra software that only necessary for a specific function by yourself.

<font style="color:rgb(52, 73, 94);">Dependencies</font> of PGAP2 are list below, and PGAP2 will check them whether in environment path or in pgap2/dependencies folder.

### Preprocessing
+ One of clustering software
    - [cd-hit](https://github.com/weizhongli/cdhit)
    - [MMseqs2](https://github.com/soedinglab/MMseqs2)
+ One of alignment software
    - [diamond](https://github.com/bbuchfink/diamond)
    - [blast+ ](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)

### Main
+ One of clustering software
    - [cd-hit](https://github.com/weizhongli/cdhit)
    - [MMseqs2](https://github.com/soedinglab/MMseqs2)
+ [mcl](https://github.com/micans/mcl)
+ One of alignment software
    - [diamond](https://github.com/bbuchfink/diamond)
    - [blast+ ](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
+ Using `--retrieve` to retrieve missing gene loci
    - [miniprot](https://github.com/lh3/miniprot)
    - [seqtk](https://github.com/lh3/seqtk)
+ Using `--reannot` to re-annotate your genome
    - [prodigal](https://github.com/hyattpd/Prodigal)

### Postprocessing
+ One of MSA software
    - [muscle](https://github.com/rcedgar/muscle)
    - [mafft](https://github.com/GSLBiotech/mafft)
    - [tcoffee](https://github.com/cbcrg/tcoffee)
+ [ClipKIT](https://github.com/JLSteenwyk/ClipKIT)
+ One of phylogenetic tree construction software
    - [IQ-TREE](http://www.iqtree.org/)
    - [FastTree](https://morgannprice.github.io/fasttree/)
    - [RAxML-ng](https://github.com/amkozlov/raxml-ng)
+ [ClonalFrameML](https://github.com/xavierdidelot/ClonalFrameML)
+ [maskrc-svg](https://github.com/kwongj/maskrc-svg)
+ [fastbaps](https://github.com/gtonkinhill/fastbaps)



### Visulization in  Preprocessing and Postprocessing modules
PGAP2 will call Rscript in your environment virable. The library should have:

+ ggpubr
+ ggrepel
+ dplyr
+ tidyr
+ patchwork
+ optparse



## Detailed documentation
Please refer documentation from [wiki](https://github.com/bucongfan/PGAP2/wiki).



