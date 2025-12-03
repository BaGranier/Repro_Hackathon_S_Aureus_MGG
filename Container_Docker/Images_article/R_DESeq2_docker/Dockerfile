# Image de base R 3.4.1
FROM rocker/r-ver:3.4.1

# Mettre à jour les dépôts Debian obsolètes et installer les dépendances système
RUN sed -i 's/deb.debian.org/archive.debian.org/g' /etc/apt/sources.list && \
    sed -i 's/security.debian.org/archive.debian.org/g' /etc/apt/sources.list && \
    sed -i '/stretch-updates/d' /etc/apt/sources.list && \
    echo 'Acquire::Check-Valid-Until "false";' > /etc/apt/apt.conf.d/99disable-check-valid-until && \
    apt-get update && apt-get install -y --no-install-recommends \
        libcurl4-openssl-dev \
        libxml2-dev \
        libssl-dev \
        libblas-dev \
        liblapack-dev \
        gfortran \
        libc6-dev \
        libncurses5-dev \
        libreadline-dev \
        libxt-dev \
        zlib1g-dev \
        libbz2-dev \
        ca-certificates \
        tar \
        bzip2 \
        wget \
    && rm -rf /var/lib/apt/lists/*

# Copier le dossier packages contenant tous les .tar.gz
COPY packages /tmp/packages

# Installer les packages R dans l'ordre correct
RUN R CMD INSTALL /tmp/packages/bitops_1.0-6.tar.gz && \
    R CMD INSTALL /tmp/packages/digest_0.6.12.tar.gz && \
    R CMD INSTALL /tmp/packages/Rcpp_0.12.12.tar.gz && \
    R CMD INSTALL /tmp/packages/RcppArmadillo_0.7.900.2.0.tar.gz && \
    \
    R CMD INSTALL /tmp/packages/plyr_1.8.4.tar.gz && \
    R CMD INSTALL /tmp/packages/reshape2_1.4.3.tar.gz && \
    R CMD INSTALL /tmp/packages/data.table_1.10.4-1.tar.gz && \
    \
    # Bloc graphique : ggplot2 et dépendances
    R CMD INSTALL /tmp/packages/colorspace_1.3-2.tar.gz && \
    R CMD INSTALL /tmp/packages/munsell_0.4.3.tar.gz && \
    R CMD INSTALL /tmp/packages/RColorBrewer_1.1-2.tar.gz && \
    R CMD INSTALL /tmp/packages/dichromat_2.0-0.tar.gz && \
    R CMD INSTALL /tmp/packages/labeling_0.3.tar.gz && \
    R CMD INSTALL /tmp/packages/scales_0.4.1.tar.gz && \
    R CMD INSTALL /tmp/packages/rlang_0.1.tar.gz && \
    R CMD INSTALL /tmp/packages/tibble_1.3.4.tar.gz && \
    R CMD INSTALL /tmp/packages/lazyeval_0.2.1.tar.gz && \
    R CMD INSTALL /tmp/packages/Formula_1.2-1.tar.gz && \
    R CMD INSTALL /tmp/packages/gtable_0.2.0.tar.gz && \
    R CMD INSTALL /tmp/packages/gridExtra_2.2.1.tar.gz && \
    R CMD INSTALL /tmp/packages/ggplot2_2.2.1.tar.gz && \
    R CMD INSTALL /tmp/packages/viridisLite_0.2.0.tar.gz && \
    R CMD INSTALL /tmp/packages/viridis_0.4.0.tar.gz && \
    \
    # Dépendances pour knitr
    R CMD INSTALL /tmp/packages/evaluate_0.10.1.tar.gz && \
    R CMD INSTALL /tmp/packages/highr_0.6.tar.gz && \
    R CMD INSTALL /tmp/packages/mime_0.5.tar.gz && \
    R CMD INSTALL /tmp/packages/markdown_0.8.tar.gz && \
    R CMD INSTALL /tmp/packages/yaml_2.1.15.tar.gz && \
    \
    # knitr + dépendances htmlTable
    R CMD INSTALL /tmp/packages/backports_1.1.0.tar.gz && \
    R CMD INSTALL /tmp/packages/knitr_1.16.tar.gz && \
    R CMD INSTALL /tmp/packages/jsonlite_1.4.tar.gz && \
    R CMD INSTALL /tmp/packages/htmltools_0.3.6.tar.gz && \
    R CMD INSTALL /tmp/packages/checkmate_1.8.3.tar.gz && \
    R CMD INSTALL /tmp/packages/htmlwidgets_0.9.tar.gz && \
    \
    # htmlTable et dépendances générales
    R CMD INSTALL /tmp/packages/htmlTable_1.9.tar.gz && \
    R CMD INSTALL /tmp/packages/base64enc_0.1-2.tar.gz && \
    R CMD INSTALL /tmp/packages/RCurl_1.95-4.10.tar.gz && \
    R CMD INSTALL /tmp/packages/latticeExtra_0.6-28.tar.gz && \
    R CMD INSTALL /tmp/packages/acepack_1.3-3.3.tar.gz && \
    \
    # Dépendances pour RSQLite
    R CMD INSTALL /tmp/packages/bit_1.1-12.tar.gz && \
    R CMD INSTALL /tmp/packages/bit64_0.9-7.tar.gz && \
    R CMD INSTALL /tmp/packages/blob_1.1.0.tar.gz && \
    R CMD INSTALL /tmp/packages/memoise_1.1.0.tar.gz && \
    R CMD INSTALL /tmp/packages/plogr_0.1-1.tar.gz && \
    R CMD INSTALL /tmp/packages/pkgconfig_2.0.1.tar.gz && \
    R CMD INSTALL /tmp/packages/BH_1.65.0-1.tar.gz && \
    R CMD INSTALL /tmp/packages/DBI_0.7.tar.gz && \
    R CMD INSTALL /tmp/packages/RSQLite_2.0.tar.gz && \
    R CMD INSTALL /tmp/packages/xtable_1.8-2.tar.gz && \
    R CMD INSTALL /tmp/packages/XML_3.98-1.9.tar.gz && \
    \
    # Bloc Bioconductor pour DESeq2
    R CMD INSTALL /tmp/packages/Hmisc_4.0-3.tar.gz && \
    R CMD INSTALL /tmp/packages/BiocGenerics_0.22.1.tar.gz && \
    R CMD INSTALL /tmp/packages/zlibbioc_1.22.0.tar.gz && \
    R CMD INSTALL /tmp/packages/S4Vectors_0.14.7.tar.gz && \
    R CMD INSTALL /tmp/packages/IRanges_2.10.5.tar.gz && \
    R CMD INSTALL /tmp/packages/XVector_0.16.0.tar.gz && \
    R CMD INSTALL /tmp/packages/GenomeInfoDbData_1.0.0.tar.gz && \
    R CMD INSTALL /tmp/packages/GenomeInfoDb_1.12.3.tar.gz && \
    R CMD INSTALL /tmp/packages/GenomicRanges_1.28.6.tar.gz && \
    R CMD INSTALL /tmp/packages/lambda.r_1.1.9.tar.gz && \
    R CMD INSTALL /tmp/packages/futile.options_1.0.0.tar.gz && \
    R CMD INSTALL /tmp/packages/futile.logger_1.4.tar.gz && \
    R CMD INSTALL /tmp/packages/snow_0.4-2.tar.gz && \
    R CMD INSTALL /tmp/packages/BiocParallel_1.10.1.tar.gz && \
    \
    R CMD INSTALL /tmp/packages/Biobase_2.36.2.tar.gz && \    
    R CMD INSTALL /tmp/packages/AnnotationDbi_1.38.2.tar.gz && \
    \
    R CMD INSTALL /tmp/packages/annotate_1.54.0.tar.gz && \    
    R CMD INSTALL /tmp/packages/genefilter_1.58.1.tar.gz && \
    R CMD INSTALL /tmp/packages/geneplotter_1.54.0.tar.gz && \
    R CMD INSTALL /tmp/packages/matrixStats_0.52.2.tar.gz && \
    R CMD INSTALL /tmp/packages/DelayedArray_0.2.7.tar.gz && \
    R CMD INSTALL /tmp/packages/SummarizedExperiment_1.6.5.tar.gz && \
    R CMD INSTALL /tmp/packages/locfit_1.5-9.1.tar.gz && \
    R CMD INSTALL /tmp/packages/DESeq2_1.16.1.tar.gz

# Définir le répertoire de travail
WORKDIR /home

# Exécuter R par défaut
CMD ["R"]

