FROM rocker/r-ver:3.4.1

LABEL maintainer="pj.gouze@gmail.com"
LABEL description="R 3.4.1 avec DESeq2 v1.7 depuis archive"

WORKDIR /app

# Installer ggplot2
RUN R -e "install.packages('ggplot2', repos='https://cloud.r-project.org')"

# Installer DESeq2 v1.7 depuis l'archive Bioconductor
RUN R -e "install.packages('https://bioconductor.org/packages/3.3/bioc/src/contrib/DESeq2_1.7.0.tar.gz', repos=NULL, type='source')"

CMD ["R"]

