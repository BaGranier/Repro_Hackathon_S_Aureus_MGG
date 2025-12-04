FROM bioconductor/bioconductor_docker:RELEASE_3_19

LABEL maintainer="vmichelet"
ENV DEBIAN_FRONTEND=noninteractive

# Paquets système nécessaires
RUN apt-get update && apt-get install -y --no-install-recommends \
        wget \
        tar \
        bzip2 \
        ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Installer BiocManager si nécessaire et DESeq2 + dépendances
RUN R -e "if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager', repos='https://cloud.r-project.org')" && \
    R -e "BiocManager::install(c('DESeq2','SummarizedExperiment','GenomicRanges','BiocParallel'))"

# Installer les packages CRAN utiles
RUN R -e "install.packages(c('data.table','ggplot2','plyr','reshape2','htmlTable','knitr','rmarkdown','jsonlite','httr','ggrepel'), repos='https://cloud.r-project.org')"

WORKDIR /home
CMD ["R"]

