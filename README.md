Putting together our model code for reproducing results of our spatiotemporal model for Malaria incidence rate prediction.

Detailed info on our model can be found in the [model directory](./model/README.md), along with notes on the dataset format and sources in the [data directory](./data/README.md).

### Reproducing

The model code was tested on an Ubuntu 22.04 container with Docker,
and a [Dockerfile](./Dockerfile) is provided for reproducing the model
environment. Typically, to install the
[CARBayesST](https://cran.r-project.org/web/packages/CARBayesST/CARBayesST.pdf)
package requires installing the build dependencies for the spatial
stack, and a verified list of packages for Ubuntu 22.04 is provided
below.

#### Installing System Packages

```sh
apt-get update
apt-get install -y --no-install-recommends \
    ca-certificates gnupg software-properties-common dirmngr wget \
    curl locales tzdata build-essential gfortran libblas-dev \
    liblapack-dev libcurl4-openssl-dev libssl-dev libxml2-dev \
    libgdal-dev libgeos-dev libproj-dev libudunits2-dev cmake
```

#### Installing R Packages

The required R packages to run the main model script
[model_analysis_batch.R](./model/model_analysis_batch.R) can be
installed with the following command:

```R
options(repos='https://cloud.r-project.org'); install.packages(c( \
    'dplyr', 'reshape2', 'ggplot2', 'MASS', 'tidyr', 'CARBayesST' \
))
```
