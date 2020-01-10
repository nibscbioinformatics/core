# Use Singularity with NextFlow

We have taken as case example the influenza workflow, because it requires a series of dependencies and in particular:

- BLAST in order to make database and search
- SEQKIT in order to convert FASTQ to FASTA
- R and numerous R packages in order to summarise the data
- Pandoc in oder to render the HTML markdown


## Creating the Singularity Image

We have created a singularity image using the following recipe, available in this folder:

```
Bootstrap: yum
OSVersion: 7
MirrorURL: http://mirror.centos.org/centos-7/7.7.1908/os/x86_64/
Include: yum

%post
yum install -y epel-release
yum install -y R
yum install -y \
    make \
    gcc \
    gcc-c++ \
    libcurl-devel \
    libxml2-devel \
    java-1.7.0-openjdk-devel \
    openssl-devel \
    texlive-* \
    pandoc
Rscript -e "install.packages(c('tidyverse','pander','rmarkdown', 'knitr'), repos = 'https://cloud.r-project.org')"
yum install -y wget
wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/ncbi-blast-2.10.0+-4.x86_64.rpm
yum install -y ncbi-blast-2.10.0+-4.x86_64.rpm
wget https://github.com/shenwei356/seqkit/releases/download/v0.11.0/seqkit_linux_amd64.tar.gz
tar -zxvf seqkit_linux_amd64.tar.gz
cp seqkit /usr/bin/.
```

**Remember you need sudo rights in order to build or edit a Singularity image**.
The image is built using the following command

```
sudo singularity build --sandbox my_image/ ~/CODE/core/singularity/basic_image.def
```

which can be further tailored if needed using the command

```
sudo singularity shell --writable my_image/
```

and then can be converted into an executable compressed image for later quick use with the following command

```
sudo singularity build my_image.sif my_image/
```

## Using the image with Nextflow


Once the image has been created (it might take a while to install everything), it can be used straight away in nextflow using the option ```-with-singularity```. However, due to the needs in mounting the appropriate use folders, it is better to use a config file, which can provide all necessary settings to make sure the working directories are visible when executing the container.
The config is also available in this folder, and is made up by the following instructions

```
process.executor = 'slurm'
process.container = '/my/path/tests/singularity/my_image.sif'
singularity.enabled = true
singularity.autoMounts = true
```
and then launched as usual, with no change in the command line except pointing to the right config file.

```
module load NextFlow/latest

nextflow \
-c ~/CODE/core/singularity/nf-singularity.config \
run ~/CODE/core/singularity/example_workflow_singularity.nf \
-with-trace \
-with-timeline nf_influenza_timeline.htm \
-with-report nf_influenza_report.htm \
--reads /your/project/raw_data \
--origin /your/project/parental.fasta \
--db "parental_db" \
--output_dir /your/path/tests/singularity/nf_test
```
