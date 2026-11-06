#!/bin/bash
# Updates the QuIPI server from the latest code on GitHub. Run on the server with:
#   bash ~/deploy_quipi/update_quipi.sh
#
# What it updates:
#   ~/deploy_quipi/   app code, quipi_env.yml, Dockerfiles, .dockerignore and this script
#   ~/run_quipi/      application.yml (the previous one is kept as application.yml.bak)
#   Docker images     quipi_base, quipi and quipi_humu
#
# The data folders are not touched and are not part of the images. ShinyProxy mounts them
# into the containers (container-volumes in application.yml):
#   ~/deploy_quipi/quipi_core/quipi_data
#   ~/deploy_quipi/quipi_humu/quipi_humu_data

# Everything runs inside main() so bash reads the whole script before starting.
# That way the script can safely replace itself with the newest version from GitHub.
main() {
    set -e

    DEPLOY=~/deploy_quipi
    RUN=~/run_quipi
    REPO=$DEPLOY/QuIPI

    # Get the latest code
    sudo rm -rf $REPO
    git clone --depth 1 https://github.com/HarrisonWismer/QuIPI.git $REPO

    # App code and environment
    cp -r $REPO/quipi_core $REPO/quipi_humu $REPO/quipi_env.yml $DEPLOY/

    # Docker build files and this script
    cp $REPO/deploy_quipi/dockerfiles/quipi_base_docker \
       $REPO/deploy_quipi/dockerfiles/quipi_docker \
       $REPO/deploy_quipi/dockerfiles/humu_docker \
       $REPO/deploy_quipi/dockerfiles/.dockerignore \
       $REPO/deploy_quipi/update_quipi.sh \
       $DEPLOY/

    # ShinyProxy config and favicon
    mkdir -p $RUN/icons
    if [ -f $RUN/application.yml ]; then
        cp $RUN/application.yml $RUN/application.yml.bak
    fi
    cp $REPO/deploy_quipi/shinyproxy/application.yml $RUN/application.yml
    cp $REPO/quipi_core/www/quipi.ico $RUN/icons/quipi.ico

    rm -rf $REPO

    # Build the images. The base image is reused from Docker's cache unless quipi_env.yml
    # changed, so a code-only update takes seconds.
    cd $DEPLOY
    sudo docker build -f quipi_base_docker -t quipi_base .
    sudo docker build -f quipi_docker -t quipi .
    sudo docker build -f humu_docker -t quipi_humu .

    # Remove the previous, now untagged, images
    sudo docker image prune -f

    # Restart ShinyProxy if it runs as a system service; otherwise ask for a manual restart
    if systemctl is-active --quiet shinyproxy; then
        sudo systemctl restart shinyproxy
        echo "Update complete. ShinyProxy restarted."
    else
        echo "Update complete. Restart ShinyProxy so it uses the new images and application.yml."
    fi
}

main "$@"
exit
