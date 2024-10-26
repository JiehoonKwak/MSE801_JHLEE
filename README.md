# 2024 MSE801 - Demo for WGS
This is only for demonstration purpose, for 2024 MSE801 KAIST  
- Reference : [GATK4 Best Practices](https://gatk.broadinstitute.org/hc/en-us), [Biostar Handbook](https://www.biostarhandbook.com/)  

We are going to use subset of data for demonstration purpose from this paper : [Human glioblastoma arises from subventricular zone cells with low-level driver mutations](https://www.nature.com/articles/s41586-018-0389-3)  

For those who wants to follow the whole workflow with real-world dataset, refer to [this data](docs/whole_setup.md) and [GATK4 Best Practices](https://gatk.broadinstitute.org/hc/en-us)


## Before we begin!
- download the data from the link provided in [HERE](docs/demo_setup.md)
```bash
bash <(curl -s https://raw.githubusercontent.com/JiehoonKwak/MSE801_JHLEE/main/download_demo.sh)
```
- install packages from `requirements.txt` in your own, new conda environment
```bash
curl -L -O https://raw.githubusercontent.com/JiehoonKwak/MSE801_JHLEE/main/requirements.txt
conda create -n YOUR_ENV -y
conda activate YOUR_ENV
conda install --file requirements.txt  -c conda-forge -c bioconda -y
```

Then, make a link to shared files from the server
```bash
# first go to your working directory
mkdir -p ref && cd ref && ln -s /home/users/SHARE/jhlee/ref/* .
```
  
Ready? Let's start!  


## Table of Contents
- Part1 : [Process Analyze-ready bam file](docs/1_pp.md)
- Part2 : [Variant Calling & Downstream Analysis](docs/2_vc.md)
- Part3 : Amplicon sequencing using `dada2` or `CRISPResso`

