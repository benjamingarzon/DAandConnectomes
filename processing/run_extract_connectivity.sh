#!/bin/sh
#source activate DAandConnectomes
#conda install -c conda-forge nilearn 
#conda install scikit-learn

WD=/home/benjamin.garzon/Data/DAD/processed/fmriprepres/fmriprep
ATLAS=/home/benjamin.garzon/Data/DAD/parcellations/shen/fconn_atlas_150_2mm.nii


OUTDIR=/home/benjamin.garzon/Data/DAD/processed/fmriprepres/connectomes/

mkdir $OUTDIR

for TASK in GNG TAB RS; do

  echo ${TASK}
  mkdir $OUTDIR/${TASK}

  for SUB in $WD/sub-*.html; do

    SUBJECT=`basename $SUB .html`
    echo $SUBJECT
    mkdir $OUTDIR/${TASK}/$SUBJECT

    INPUT=$WD/${SUBJECT}/func/${SUBJECT}_task-${TASK}rest_run-1_bold_space-MNI152NLin2009cAsym_variant-smoothAROMAnonaggr_preproc.nii.gz
    OUTPUT=$OUTDIR/${TASK}/$SUBJECT/zFC_150.csv
    CONFOUNDS=$WD/${SUBJECT}/func/${SUBJECT}_task-${TASK}rest_run-1_bold_confounds.tsv
    GLOBAL=$WD/${SUBJECT}/func/${SUBJECT}_task-${TASK}rest_GlobalSignal.tsv
    # extract only relevant columns
#    awk '{print $3}' $CONFOUNDS | tail -n +2 > $GLOBAL
    awk '{print $3}' $CONFOUNDS > $GLOBAL

    ./extract_connectivity.py $INPUT $ATLAS $GLOBAL $OUTPUT 2.0 full 0
  done
done
