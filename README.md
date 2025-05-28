## Introduction
BioMGCore is designed to determine the biological metabolic gene clusters in microorganisms. The input files are expected to involve the bacterial genomes, or the species-level genome bins (SGBs). Currently not applicable to fungi and other higher organisms.

## Dependencies
gutSMASH: https://gutsmash.bioinformatics.nl/
antiSMASH: https://antismash.secondarymetabolites.org/
Python3: https://www.python.org/downloads/

### Repository
#### Step 1. Install the Dependencies
The [Dependencies](#Dependencies) are required to be installed and added to the system `$PATH`
#### Step 2. Download the BioMGCore toolkit here

## Running
### 1. Genome annotation by gutSMASH and antiSMASH
#### 1-1. gutSMASH annotation:
```
python3  run_gutsmash.py --cb-knownclusters --enable-genefunctions  your_data_path/genome_file  --output-dir  examples/caiTABCDE/results_genome_file_name  --genefinding-tool  prodigal  -c  16
```
#### 1-2. antiSMASH annotation:
```
antismash --taxon bacteria --output-dir examples/thiopeptide/results_genome_file_name --databases your_data_path/antismash/databases --genefinding-tool prodigal-m --cb-general --cb-knownclusters --cb-subclusters --asf --pfam2go --smcog-trees --cpus 10 your_data_path/genome_file
```

### 2. Sequences extraction
#### 2-1. Sequence extraction for gutSMASH (example):
```
python3  BioMGCore.py  -n  caiTABCDE  -i  examples/caiTABCDE  -o output_caiTABCDE -g 5 -t s
```
#### 2-2. Sequence extraction for antiSMASH (example):
```
python3  BioMGCore.py  -n  thiopeptide  -i  examples/thiopeptide  -o  output_thiopeptid
```

### 3. Statistics for antiSMASH annotations:
```
python3  antiSTAT.py -i  data_path_antiSMASH_annotation -o  ./statistics.xlsx
```


## Copyrights
Shulei Jia: jiashulei@tmu.edu.cn

School of Basic Medical Sciences, Tianjin Medical University, Tianjin, 300070, China

Guang Yang: 2020000050@jou.edu.cn

Jiangsu Marine Resources Development Research Institute, Jiangsu Ocean University, Lianyungang, 222005, China
