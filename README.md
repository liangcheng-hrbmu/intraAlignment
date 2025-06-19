# Requirements

### softwares and packages requirements

```bash
name: intraAlignment
channels:
  - bioconda
  - conda-forge
  - defaults
dependencies:
  - blast
  - bracken=2.9=py310h0dbaff4_0
  - kraken2=2.1.3=pl5321hdcf5f25_0
  - python=3.10.11=h7a1cb2a_2
  - samtools
  - taxonkit
  - pandas
  - tqdm
  - colorama
```
# Installation
```bash
conda env create -n intraAlignment -f intraAlignment.yaml
```



# Example

```bash
python [...]/intraAlignment.py -rp [...]/cellranger[&spaceranger]/outs -method kraken2 -db [...]/k2_db


python [...]/intraAlignment.py -rp [...]/cellranger[&spaceranger]/outs -method blast -db [...]/nt_db/nt
```








