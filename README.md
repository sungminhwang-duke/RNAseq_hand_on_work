# RNAseq_hand_on_work

Part I - for processing the raw sequencing data <br/>
1. Go to the folder "Seup_env" and download a file to make a virtual environment for a data process. <br/>
2. Execute this command in the linux-based operational system: <br/>
    For Mac: _conda env create -f condaTESTmac.yml_<br/>
    For Linux: _conda env create -f condaTESTlinux.yml_<br/>
   -Keep in mind that the downloaded file path should be identical to the Terminal, "condaTESTmac.yml". <br/>
   -Activate the env with this command: _source activate condaTESTmac_ <br/>
3. Download raw data (toy data) and materials, as well as run the code with this command "_bash run-rnaseq.sh_"  <br/>

Part II - for the data analysis
1. For the DEG analysis, go to "DEGs" folder and download the exercise dataset (fastq.gz of toy data) and DESeq2 command lines. <br/>
2. Execute the DESeq2 command in R (Windows, Mac, etc...)
