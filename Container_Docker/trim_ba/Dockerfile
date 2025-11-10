# syntax=docker/dockerfile:1
FROM --platform=linux/amd64 ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive
WORKDIR /opt

# 1️⃣ Install deps (wget, bzip2, gcc, make…)
RUN apt-get update && apt-get install -y \
    wget bzip2 build-essential ca-certificates && \
    rm -rf /var/lib/apt/lists/*

# 2️⃣ Install Miniconda2 manually
RUN wget https://repo.anaconda.com/miniconda/Miniconda2-latest-Linux-x86_64.sh -O miniconda.sh && \
    bash miniconda.sh -b -p /opt/conda && rm miniconda.sh
ENV PATH=/opt/conda/bin:$PATH

# 3️⃣ Create conda env + install Cutadapt 1.11 + Trim Galore 0.6.6
RUN conda create -y -n cutadapt_1_11 python=2.7 "cython<0.29" setuptools wget && \
    bash -c "source activate cutadapt_1_11 && \
      wget https://github.com/marcelm/cutadapt/archive/refs/tags/v1.11.tar.gz && \
      tar -xzf v1.11.tar.gz && cd cutadapt-1.11 && python setup.py install && \
      cd /opt && rm -rf v1.11.tar.gz cutadapt-1.11 && \
      wget https://github.com/FelixKrueger/TrimGalore/archive/refs/tags/0.6.6.tar.gz -O trim_galore-0.6.6.tar.gz && \
      tar -xzf trim_galore-0.6.6.tar.gz && cd TrimGalore-0.6.6 && chmod +x trim_galore && \
      ln -s $(pwd)/trim_galore /opt/conda/envs/cutadapt_1_11/bin/trim_galore && \
      cd /opt && rm -rf trim_galore-0.6.6.tar.gz"

# 4️⃣ Environment par défaut
SHELL ["/bin/bash", "-c"]
ENV PATH=/opt/conda/envs/cutadapt_1_11/bin:$PATH
WORKDIR /data

# 5️⃣ (optionnel mais utile) Vérifie Trim Galore à la fin du build
RUN trim_galore --version || true

ENTRYPOINT ["trim_galore"]
