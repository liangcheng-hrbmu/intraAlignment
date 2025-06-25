# 1. Requirements

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

# 2. Installation

### **method 1 - git & conda**

```bash
git clone https://github.com/liangcheng-hrbmu/intraAlignment.git
cd intraAlignment
conda env create -n intraAlignment -f intraAlignment.yaml
```

### **method 2 - docker**

```bash
docker pull ghcr.io/liangcheng-hrbmu/intra_alignment:latest
```



# 3. Reference database preparation

**We recommend selecting a disk location with sufficient available space to download and store the reference databases.**

1. **[kraken2 reference database](https://benlangmead.github.io/aws-indexes/k2)**

****

```bash
cd <ref_db_path> # disk location with sufficient available space
wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20250402.tar.gz # 86.8 GB
mkdir -p k2_standard_20250402 && tar -zxvf k2_standard_20250402.tar.gz -C k2_standard_20250402 #
```



2. **[blastn reference database](https://ftp.ncbi.nlm.nih.gov/blast/db/)**

****

```bash
# option 1
cd <ref_db_path> # disk location with sufficient available space
update_blastdb.pl --passive --decompress nt

# option 2
# curl -s ftp://ftp.ncbi.nlm.nih.gov/blast/db/ | grep -o 'nt\.[0-9]\{3\}\.tar\.gz\(\.md5\)\?' | awk '{print "https://ftp.ncbi.nlm.nih.gov/blast/db/" $0}' - > file_url.txt
# cat -s file_list.txt | grep -o 'nt\.[0-9]\{3\}\.tar\.gz\(\.md5\)\?' | awk '{print "https://ftp.ncbi.nlm.nih.gov/blast/db/" $0}' - > file_url.txt
# aria2c -i file_url.txt --enable-http-pipelining="true" -x16 -s100 -j4 --allow-piece-length-change="true" --piece-length="16M" --min-split-size="16M" -c -m 5 --retry-wait=20
```

   



# 4. Usage example

### **method 1 - git & conda**

**10X data**

```bash
# 10X using kraken2 
python <...>/intraAlignment.py -rp <10x_dt_path>/<cellranger[&spaceranger]>/outs -db <ref_db_path>/k2_standard_20250402

# 10X usingblastn 
python <...>/intraAlignment.py -rp <10x_dt_path>/<cellranger[&spaceranger]>/outs -method blast -db <ref_db_path>/nt_db/nt
```

**C4 data**

```bash
# C4 usingkraken2 
python [...]/intraAlignment.py -p c4 -c4_o <c4_dt_path>/<dnbc4tools>/outs -c4_r1 <c4_dt_path>/<sample>_cDNA_R1.fq.gz -db <ref_db_path>/k2_standard_20250402

# C4 using blastn 
python [...]/intraAlignment.py -p c4 -c4_o <c4_dt_path>/<dnbc4tools>/outs -c4_r1 <c4_dt_path>/<sample>_cDNA_R1.fq.gz -method blast -db <ref_db_path>/nt_db/nt
```



### **method 2 - docker**

```bash
# create contanier
docker run -it --name intraAlignment -v <data_path>:/<data_path> ghcr.io/liangcheng-hrbmu/intra_alignment

intraAlignment -rp <10x_dt_path>/<cellranger[&spaceranger]>/outs -db <ref_db_path>/k2_standard_20250402
```

