# Use Ubuntu (Jammy 22.04 is stable with good package availability)
FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

# 1) Base OS deps and tools
RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates \
    gnupg \
    software-properties-common \
    dirmngr \
    wget \
    curl \
    locales \
    tzdata \
    build-essential \
    gfortran \
    libblas-dev \
    liblapack-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    # spatial stack commonly required by sp/sf/spdep/CARBayes*:
    libgdal-dev \
    libgeos-dev \
    libproj-dev \
    libudunits2-dev \
    cmake

# Set a UTF-8 locale (helps some R packages / plotting)
RUN locale-gen en_US.UTF-8
ENV LANG=en_US.UTF-8
ENV LC_ALL=en_US.UTF-8

# 2) Add the official CRAN repo for R (for Ubuntu 22.04 "jammy")
# See: https://cran.r-project.org/bin/linux/ubuntu/
RUN mkdir -p /etc/apt/keyrings && \
    curl -fsSL https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc \
      | gpg --dearmor -o /etc/apt/keyrings/cran-archive-keyring.gpg && \
    echo "deb [signed-by=/etc/apt/keyrings/cran-archive-keyring.gpg] https://cloud.r-project.org/bin/linux/ubuntu jammy-cran40/" \
      > /etc/apt/sources.list.d/cran-r.list

# 3) Install R
RUN apt-get update && apt-get install -y --no-install-recommends \
    r-base r-base-dev \
    && rm -rf /var/lib/apt/lists/*

# 4) Preinstall required R packages
# Note: 'splines' ships with base R; no need to install it from CRAN.
# We explicitly install the rest from a reliable CRAN mirror.
RUN R -q -e "options(repos='https://cloud.r-project.org'); install.packages(c( \
    'dplyr', 'reshape2', 'ggplot2', 'MASS', 'tidyr', 'CARBayesST' \
))"

# 5) App directory and script
RUN R -q -e "options(repos='https://cloud.r-project.org'); install.packages(c('CARBayesST'), Ncpus = 4)"

COPY model /app/model
COPY data /app/data

WORKDIR /app/model

# 6) Default command: run the script
CMD ["Rscript", "model_analysis_batch.R"]
