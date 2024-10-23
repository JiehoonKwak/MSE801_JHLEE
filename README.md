# 2024 MSE801 - Demo for WGS
This is only for demonstration purpose, for 2024 MSE801 KAIST  
- Reference : [GATK4 Best Practices](https://gatk.broadinstitute.org/hc/en-us), [Biostar Handbook](https://www.biostarhandbook.com/)  

We are going to use subset of data for demonstration purpose from this paper : [Human glioblastoma arises from subventricular zone cells with low-level driver mutations](https://www.nature.com/articles/s41586-018-0389-3)  

For those who wants to follow the whole workflow with real-world dataset, refer to [this data](docs/setup.md) and [whole workflow](docs/whole_workflow.md)


## Before we begin!
- download the data from the link provided in [HERE](docs/demo_setup.md)
- install packages from `requirements.txt` in your own, new conda environment
```bash
conda create -n YOUR_ENV -y
conda activate YOUR_ENV
conda install --file requirements.txt
```

## Table of Contents
- Part1 : [Process Analyze-ready bam file](docs/1_pp.md)
- Part2 : [Variant Calling & Downstream Analysis](docs/2_vc.md)
- Part3 : Amplicon sequencing using `dada2`

