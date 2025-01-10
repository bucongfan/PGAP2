Pan-Genome Analysis Pipeline 2



## Quick start
### [<font style="color:rgb(51, 51, 51);">Basic usage</font>](https://gthlab.au/panaroo/#/gettingstarted/quickstart?id=basic-usage)
The inputdir contains all genome and its annotation file. 

PGAP2 supports multiple format as input, <font style="color:rgb(52, 73, 94);">GFFs in the same format as output by Prokka, GFFs and its Genome fasta file in different file, Genbank flat file (GBFF) or just Genome fasta file (--reannot is required).</font>

<font style="color:rgb(52, 73, 94);">Each format of input file can mix in one input directory, PGAP2 will recognize and process them according to its prefix and suffix.</font>

```bash
pgap2 main -i inputdir/ -o outputdir/
```

### Preprocessing
Quality check and visulization conducted by PGAP2 in preprocessing step. PGAP2 will generate a interactivae HTML file and coresponding vector figures to direct user understanding its input data. Input data and pre-alignment result will stroe as a pickle file for quick re-start the same calculation step.

```bash
pgap2 prep -i inputdir/ -o outputdir/
```

### Postprocessing
Postprocessing pipeline performed by PGAP2. There are various pipeline integrated in postprocessing module, such as statistical analysis, single copy tree building, population clustering and tajima's D test. Whatever which submodule you want to carry out, you can always run as:

```bash
pgap2 post [submodule] [options] -i inputdir/ -o outputdir/
```

The inputdir is the outputdir of main module.

PGAP2 also support statistical analysis using a PAV file indepandently:

```bash
pgap2 post profle --pav your_pav_file -o outputdir/
```

## Installation
The best way to install full version of PGAP2 package is using conda:

```bash
conda create -n pgap2 -c conda-forge -c bioconda -c defaults pgap2
```

<font style="color:rgb(52, 73, 94);">alternatively it is often faster to use the </font>[mamba](https://github.com/mamba-org/mamba)<font style="color:rgb(52, 73, 94);"> solver:</font>

```bash
conda create -n pgap2 -c conda-forge mamba
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
+ [cd-hit](about:blank)
+ One of alignment software
    - [diamond](https://github.com/bbuchfink/diamond)
    - [blast+ ](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)

### Main
+ [cd-hit](about:blank)
+ [mcl](https://github.com/micans/mcl)
+ One of alignment software
    - [diamond](https://github.com/bbuchfink/diamond)
    - [blast+ ](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
+ If need to retrive missing gene locus with --retrieve
    - [miniprot](https://github.com/lh3/miniprot)
    - [seqtk](https://github.com/lh3/seqtk)
+ If need to re-annot your genome with --reannot
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



## Detilled <font style="color:rgb(31, 35, 40);">Documentation</font>
Please refer detailled documentation from here
