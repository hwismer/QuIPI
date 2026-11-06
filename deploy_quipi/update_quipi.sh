#!/bin/bash
# Pulls the latest QuIPI code from GitHub and rebuilds the Docker images.
# Run on the server with:  bash ~/deploy_quipi/update_quipi.sh
#
# The data folders are not part of the images. ShinyProxy mounts them into the containers
# (container-volumes in application.yml), so they stay where they are:
#   ~/deploy_quipi/quipi_core/quipi_data
#   ~/deploy_quipi/quipi_humu/quipi_humu_data
set -e

DEPLOY=~/deploy_quipi
cd $DEPLOY

# Copy the latest code and deploy files from GitHub
sudo rm -rf QuIPI
git clone --depth 1 https://github.com/HarrisonWismer/QuIPI.git
cp -r QuIPI/quipi_core QuIPI/quipi_humu QuIPI/quipi_env.yml $DEPLOY/
cp QuIPI/deploy_quipi/dockerfiles/quipi_base_docker QuIPI/deploy_quipi/dockerfiles/quipi_docker \
   QuIPI/deploy_quipi/dockerfiles/humu_docker QuIPI/deploy_quipi/dockerfiles/.dockerignore $DEPLOY/
rm -rf QuIPI

# The base image is only rebuilt when quipi_env.yml changes; the app images take seconds
sudo docker build -f quipi_base_docker -t quipi_base .
sudo docker build -f quipi_docker -t quipi .
sudo docker build -f humu_docker -t quipi_humu .

# Remove the previous, now untagged, images
sudo docker image prune -f

echo "Images rebuilt. Restart ShinyProxy so it starts containers from the new images."
