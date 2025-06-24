#!/bin/bash

# Update commands
sudo apt update
sudo apt-get update -y
sudo apt-get dist-upgrade -y
# Install some system libraries
sudo apt-get -y install \
    nginx \
    gdebi-core \
    apache2-utils \
    pandoc \
    pandoc-citeproc \
    libssl-dev \
    libcurl4-gnutls-dev \
    libcairo2-dev \
    libgsl0-dev \
    libgdal-dev \
    libgeos-dev \
    libproj-dev \
    libxml2-dev \
    libxt-dev \
    libv8-dev \
    libhdf5-dev \
    git \
    gh

sudo apt install build-essential

sudo apt update -qq

sudo apt install --no-install-recommends software-properties-common dirmngr

wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | sudo tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc

sudo add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"

sudo apt-get install r-base

sudo su - -c "R -e \"install.packages('shiny', repos='https://cran.rstudio.com/')\""

sudo apt-get install gdebi-core

wget https://download3.rstudio.org/ubuntu-18.04/x86_64/shiny-server-1.5.22.1017-amd64.deb

sudo gdebi shiny-server-1.5.22.1017-amd64.deb
