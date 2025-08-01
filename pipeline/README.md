# Pipeline README
-----------------------------------------------------
By Alexander Balzer, August 2025
This project is part of my bachelor thesis and is used to to visualize the HGT score of amino acids of the N-terminal end 
of all proteins belonging to the mitochondria having a cleavable mitochondrial targeting signal.
-----------------------------------------------------
## Steps
Proteom(fasta) + Annotations (goa)
| filtered by GO term
V
mitochondrial proteins
| filter by cleavable MTS using MitoFates as an api
V
mitochondrial proteins with
cleavable MTS
[HGT/absolute values]    -> visualize as heatmap

## Usage
1. choose a input folder and provide it with:
    - GO annotations from https://ftp.ebi.ac.uk/pub/databases/GO/goa/proteomes/ 
    - Proteomes from https://www.uniprot.org/proteomes/
    Both files have to have the same name, the GO annotations with .goa and the Proteome with .fasta
    the Names should look something like saccharomyces_cerevisiae.fasta and saccharomyces_cerevisiae.goa
    as they are used later as labels in the heatmap.
2. configure the config file config.ini and provide it with your input directory and the link to the MitoFates perl script.
3. write into the flaglist.txt file the name of your organism and if its a plant, a fungi or a metazoa
    use the same name as under point 1.
    it should look like this:
    saccharomyces_cerevisiae: fungi
4. Run the pipeline with:
    ```bash
    python main.py
    ```

## Directory Structure
```
pipeline/
├── input/           # Input data
├── output/          # output data
├── cache/           # cache data
├── config.ini       # configuration file
├── README.md        # Documentation
├── flaglist.txt     # flags for MitoFates
└── bachelor_env.yml # conda environment for the project
```

Aditional Information:
If you use MitoFates for your research, please cite their paper.
Website: https://mitf.cbrc.pj.aist.go.jp/MitoFates/cgi-bin/top.cgi

Paper:
MitoFates: Improved Prediction of Mitochondrial Targeting Sequences and Their Cleavage Sites.
Yoshinori Fukasawa, Junko Tsuji, Szu-Chin Fu, Kentaro Tomii, Paul Horton and Kenichiro Imai.
Molecular & Cellular Proteomics, 14(4): 1113-1126, 2015
