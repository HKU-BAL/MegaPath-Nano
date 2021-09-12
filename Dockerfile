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
    conda install pandas psutil pybedtools porechop==0.2.4 bioconvert seqtk minimap2 bcftools samtools==1.9 'pysam>=0.16.0' tabulate cgecore==1.5.6 'ncbi-amrfinderplus>=3' 'rgi>=5' clair=2.1.1 parallel=20191122"
RUN cd MegaPath-Nano/bin/realignment/realign/ && \ 
    g++ -std=c++14 -O1 -shared -fPIC -o realigner ssw_cpp.cpp ssw.c realigner.cpp && \
    g++ -std=c++11 -shared -fPIC -o debruijn_graph -O3 debruijn_graph.cpp && \
    gcc -Wall -O3 -pipe -fPIC -shared -rdynamic -o libssw.so ssw.c ssw.h && \
    cd - && \
    cd MegaPath-Nano/bin/Clair-ensemble/Clair.beta.ensemble.cpu/clair/ && \
    g++ ensemble.cpp -o ensemble
    cd - && \
    cd MegaPath-Nano/bin/samtools-1.13 && \
    ./configure && make && make install


