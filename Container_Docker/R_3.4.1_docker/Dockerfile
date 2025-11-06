FROM rocker/r-ver:3.4.1

LABEL maintainer="pj.gouze@gmail.com"
LABEL description="Image R 3.4.1 basée sur rocker/r-ver"

# Définir le répertoire de travail
WORKDIR /app

# (Optionnel) Installer des packages R
RUN R -e "install.packages('ggplot2', repos='https://cloud.r-project.org')"

# Commande par défaut
CMD ["R"]

