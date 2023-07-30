FROM gcr.io/ris-registry-shared/rstudio:4.2.3

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y libsodium-dev samba-client vsftpd libprotobuf-dev protobuf-compiler cmake libfribidi-dev libharfbuzz-dev libfontconfig1-dev libcurl4-openssl-dev libssl-dev libxml2-dev libudunits2-dev libgdal-dev 

# Get and install system dependencies
RUN R -e "install.packages('remotes', repos='http://cran.us.r-project.org')" && R -e "remotes::install_github('r-hub/sysreqs')"

ENV RENV_VERSION 0.17.3
RUN R -e "install.packages('remotes', repos = c(CRAN = 'https://cloud.r-project.org'))"
RUN R -e "remotes::install_github('rstudio/renv@${RENV_VERSION}')"

WORKDIR .
COPY renv.lock renv.lock

RUN R -e 'renv::restore()'