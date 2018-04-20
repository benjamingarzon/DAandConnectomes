#!/bin/sh

# create singularity file container
docker run --privileged -t --rm \
    -v /var/run/docker.sock:/var/run/docker.sock \
    -v /tmp/fmriprep:/output \
    singularityware/docker2singularity \
    poldracklab/fmriprep:latest

cp /tmp/fmriprep/poldracklab_fmriprep_latest-*.img /home/benjamin.garzon/Data/DAD/processed/fmriprep.img

WD=/home/benjamin.garzon/Data/DAD/processed/fmriprep
WORK=/home/benjamin.garzon/Data/DAD/processed/work
rm -r $WD
rm -r $WORK

mkdir $WD
cp /usr/local/freesurfer/license.txt $WD/

PYTHONPATH="" singularity run /home/benjamin.garzon/Data/DAD/processed/fmriprep.img /home/benjamin.garzon/Data/DAD/data_BIDS $WD participant --use-aroma --use-syn-sdc --fs-license-file $WD/license.txt --nthreads 20 -w $WORK --omp-nthreads 20
 #--participant-label --write-graph


