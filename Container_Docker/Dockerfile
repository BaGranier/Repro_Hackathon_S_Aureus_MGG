FROM ubuntu:22.04

LABEL maintainer="pierje.goosse@example.com"
LABEL description="Container Docker pour le SRA Toolkit"

RUN apt-get update && apt-get install -y \
    wget \
    tar \
    gzip \
    && rm -rf /var/lib/apt/lists/*

RUN wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz && \
    tar -xzf sratoolkit.current-ubuntu64.tar.gz && \
    mv sratoolkit.* /opt/sratoolkit && \
    ln -s /opt/sratoolkit/bin/* /usr/local/bin/ && \
    rm sratoolkit.current-ubuntu64.tar.gz

WORKDIR /data

VOLUME ["/data"]

CMD ["fastq-dump", "--version"]

