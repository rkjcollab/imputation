# Container for the pipeline, used by both Docker and apptainer. Should be
# pulled from GHCR.
# If need to re-build, run the following from the top level of this repository
# with:
#   docker buildx build --platform linux/amd64 -t <image>:<tag> - < envs/Dockerfile
FROM python:3.11-slim-bookworm

LABEL author="SDS"

# General tools, bcftools build dependencies, and R
RUN apt-get update -y && apt-get install -y \
        wget \
        unzip \
        g++ \
        make \
        bzip2 \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
        libcurl4-openssl-dev \
        r-base \
        libxml2-dev \
        vcftools \
        7zip \
        curl \
        graphviz \
        openjdk-17-jdk \
        aria2 \
        jq \
    && rm -rf /var/lib/apt/lists/*

# CrossMap, Snakemake, and Pulp (pinned to 2.7.0 due to list solvers bug)
RUN pip3 install pip==26.2.2 && \
    pip3 install CrossMap==0.7.4 requests==2.34.2 snakemake==9.27.1 pulp==2.7.0

# R packages from a dated CRAN snapshot
RUN Rscript -e 'install.packages(c("argparse"), repos="https://packagemanager.posit.co/cran/2026-11-01")'

# bcftools, and tabix and bgzip from its bundled htslib, built from source
# since the Debian version is older
RUN wget https://github.com/samtools/bcftools/releases/download/1.24/bcftools-1.24.tar.bz2 && \
    tar -xjf bcftools-1.24.tar.bz2 && \
    cd bcftools-1.24 && \
    ./configure --prefix=/usr/local && \
    make && \
    make install && \
    cd htslib-1.24 && \
    make install && \
    cd ../.. && \
    rm -r bcftools-1.24 bcftools-1.24.tar.bz2

# plink 1.9
RUN wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20260927.zip && \
    unzip plink_linux_x86_64_20260927.zip -d plink_linux_x86_64/ && \
    mv plink_linux_x86_64/plink /usr/local/bin && \
    rm plink_linux_x86_64_20260927.zip && \
    rm -r plink_linux_x86_64

# plink 2
RUN wget https://s3.amazonaws.com/plink2-assets/alpha7/plink2_linux_x86_64_20261015.zip && \
    unzip plink2_linux_x86_64_20261015.zip -d plink2_linux_x86_64/ && \
    mv plink2_linux_x86_64/plink2 /usr/local/bin && \
    rm plink2_linux_x86_64_20261015.zip && \
    rm -r plink2_linux_x86_64

# Check tool versions at build time
RUN plink --version && \
    plink2 --version && \
    CrossMap --version && \
    R --version && \
    python --version && \
    bcftools --version | head -n 1 && \
    tabix --version | head -n 1 && \
    bgzip --version | head -n 1 && \
    snakemake --version && \
    7zz | head -n 2 && \
    Rscript -e 'library(argparse)'
