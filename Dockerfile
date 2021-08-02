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
        git && \
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
    conda create -n mpn python=3.6.10
RUN echo "source activate mpn" > ~/.bashrc
ENV PATH /opt/conda/envs/mpn/bin:$PATH
RUN /bin/bash -c ". activate mpn && \
    conda install pandas==1.1.5 psutil==5.6.5 pybedtools==0.8.0 porechop==0.2.4 bioconvert==0.4.3 seqtk==1.3 minimap2==2.21 bcftools==1.9 samtools==1.9 pysam==0.16.0 tabulate==0.8.9 cgecore==1.5.6 ncbi-amrfinderplus==3.10.5 rgi==5.2.0 biopython==1.72 pyahocorasick==1.1.7 filetype==1.0.7libdeflate==1.6"
