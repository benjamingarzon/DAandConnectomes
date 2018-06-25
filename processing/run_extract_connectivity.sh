#!/bin/sh
#source activate DAandConnectomes
#conda install -c conda-forge nilearn 
#conda install scikit-learn

WD=/home/benjamin.garzon/Data/DAD/processed/fmriprep/fmriprep
ATLAS=/home/benjamin.garzon/Data/DAD/parcellations/shen/fconn_atlas_150_2mm.nii


OUTDIR=/home/benjamin.garzon/Data/DAD/processed/fmriprep/connectomes/

mkdir $OUTDIR

for TASK in GNG TAB RS; do

  echo ${TASK}
  mkdir $OUTDIR/${TASK}

  for SUB in $WD/sub-???; do
  
    LIST=""
    for RUN in 1 2 3; do
    
    SUBJECT=`basename $SUB .html`
    echo $SUBJECT
    mkdir $OUTDIR/${TASK}/$SUBJECT

    INPUT=$WD/${SUBJECT}/func/${SUBJECT}_task-${TASK}rest_run-${RUN}_bold_space-MNI152NLin2009cAsym_variant-smoothAROMAnonaggr_preproc.nii.gz
    OUTPUT=$OUTDIR/${TASK}/$SUBJECT/zFC_150-${RUN}.csv
    CONFOUNDS=$WD/${SUBJECT}/func/${SUBJECT}_task-${TASK}rest_run-${RUN}_bold_confounds.tsv
    GLOBAL=$WD/${SUBJECT}/func/${SUBJECT}_task-${TASK}rest_GlobalSignal.tsv
    if [ -f "$INPUT" ]; then
    awk '{print $3}' $CONFOUNDS > $GLOBAL
    # extract only relevant columns
#    awk '{print $3}' $CONFOUNDS | tail -n +2 > $GLOBAL
    ./extract_connectivity.py $INPUT $ATLAS $GLOBAL $OUTPUT 2.0 full 0
    LIST="$LIST $OUTPUT"
    fi

    done
    echo $LIST
    ./average_connectivity.py $OUTDIR/${TASK}/$SUBJECT/zFC_150.csv "$LIST"
  done
done
