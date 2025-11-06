FROM ubuntu:14.04

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    wget \
    make \
    unzip \
    zlib1g-dev gcc-4.4 g++-4.4\
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/
RUN wget http://sourceforge.net/projects/bowtie-bio/files/bowtie/0.12.7/bowtie-0.12.7-src.zip/download -0 bowtie-0.12.7-src.zip
RUN unzip bowtie-0.12.7-src.zip
RUN rm bowtie-0.12.7-src.zip

WORKDIR /opt/bowtie-0.12.7

RUN make CC=gcc-4.4 CXX=g++-4.4

ENV PATH="/opt/bowtie-0.12.7:$PATH"

CMD ["bowtie", "--version"]
