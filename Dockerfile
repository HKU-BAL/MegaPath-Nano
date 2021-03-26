FROM ubuntu:16.04

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8 PATH=/opt/MegaPath-Nano/bin:/opt/conda/bin:$PATH

# update ubuntu packages
RUN apt-get update --fix-missing && \
    yes|apt-get upgrade && \
    apt-get install -y \
        wget \
        bzip2 \
        make \
        gcc \
        g++ \
        git \
        vcftools && \
    rm -rf /bar/lib/apt/lists/* && \
    git clone --depth 1 https://github.com/HKU-BAL/MegaPath-Nano /opt/MegaPath-Nano

WORKDIR /opt/MegaPath-Nano
COPY . .

# install anaconda
RUN wget --quiet https://repo.anaconda.com/archive/Anaconda3-2020.07-Linux-x86_64.sh && \
    bash Anaconda3-2020.07-Linux-x86_64.sh -b -p /opt/conda && \
    rm Anaconda3-2020.07-Linux-x86_64.sh

# create conda environment
RUN conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda create -n mpn python=3.6
RUN echo "source activate mpn" > ~/.bashrc
ENV PATH /opt/conda/envs/mpn/bin:$PATH
RUN /bin/bash -c ". activate mpn && \
    conda install pandas==0.23 psutil pybedtools porechop bioconvert seqtk minimap2 bcftools samtools pysam tabulate cgecore ncbi-amrfinderplus && \
    pip install -Iv biopython==1.72 && \
    conda update -c conda-forge biopython && \
    pip install git+https://github.com/arpcard/rgi.git pyfaidx pyahocorasick seaborn"
