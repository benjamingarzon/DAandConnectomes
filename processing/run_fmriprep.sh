#!/bin/sh

# create singularity file container
#docker run --privileged -t --rm \
#    -v /var/run/docker.sock:/var/run/docker.sock \
#    -v /tmp/fmriprep:/output \
#    singularityware/docker2singularity \
#    poldracklab/fmriprep:latest

#cp /tmp/fmriprep/poldracklab_fmriprep_latest-*.img /home/benjamin.garzon/Data/DAD/processed/fmriprep.img

WD=/home/benjamin.garzon/Data/DAD/processed/fmriprep
WORK=/home/benjamin.garzon/Data/DAD/processed/work

rm -r $WD
rm -r $WORK
mkdir $WD
cp /usr/local/freesurfer/license.txt $WD/


for s in sub-D02 sub-D03 sub-D04 sub-D05 sub-D07 sub-D08 sub-D11 sub-D13 sub-D14 sub-D15 sub-D16 sub-D17 sub-D18 sub-D19 sub-D20 sub-D21 sub-D22 sub-D23 sub-D24 sub-D25 sub-D26 sub-D29 sub-D30 sub-D31 sub-D33 sub-D34 sub-D35 sub-D36 sub-D37 sub-D38 sub-D39 sub-D40 sub-D42 sub-D43 sub-D46 sub-D47 sub-D48 sub-D49 sub-D50 sub-D52 sub-D55 sub-D56 sub-D58 sub-D60 sub-D61 sub-D62 sub-D63 sub-D64 sub-D65 sub-D66 sub-D67 sub-D70 sub-D72 sub-D80 sub-D82 sub-D83 sub-D84 sub-D85 sub-D86 sub-D90;  
do
  WORK=/home/benjamin.garzon/Data/DAD/processed/work/
  PYTHONPATH="" singularity run /home/benjamin.garzon/Data/DAD/processed/fmriprep.img /home/benjamin.garzon/Data/DAD/data_BIDS $WD participant --use-aroma --use-syn-sdc --fs-license-file $WD/license.txt --nthreads 4 -w $WORK --participant-label $s & 
done

PYTHONPATH="" singularity run /home/benjamin.garzon/Data/DAD/processed/fmriprep.img /home/benjamin.garzon/Data/DAD/data_BIDS $WD participant --use-aroma --use-syn-sdc --fs-license-file $WD/license.txt --nthreads 4 -w $WORK --participant-label $s


