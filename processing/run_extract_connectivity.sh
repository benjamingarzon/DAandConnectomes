#!/bin/sh
#source activate DAandConnectomes
#conda install -c conda-forge nilearn 
#conda install scikit-learn

WD=/home/benjamin.garzon/Data/DAD/processed/fmriprep/fmriprep
ATLAS=/home/benjamin.garzon/Data/DAD/parcellations/shen/fconn_atlas_150_2mm.nii

MEAN_THR=0.3
MAX_THR=5
PROP_THR=0.3
MAX_INVALID=20

OUTDIR=/home/benjamin.garzon/Data/DAD/processed/fmriprep/connectomes/

mkdir $OUTDIR
echo "SUBJECT;TASK;RUN;FD_MEAN;FD_MAX;FD_2;FD_3;INVALID" > $OUTDIR/FramewiseDisplacement.csv

for TASK in GNG TAB RS; do
  rm $OUTDIR/Subjects_${TASK}.txt
  echo ${TASK}
  mkdir $OUTDIR/${TASK}

  for SUB in $WD/sub-???; do
  
    LIST=""
    for RUN in 1 2 3; do
    
    SUBJECT=`basename $SUB .html`
    echo $SUBJECT
    mkdir $OUTDIR/${TASK}/$SUBJECT

    INPUT=$WD/${SUBJECT}/func/${SUBJECT}_task-${TASK}rest_run-${RUN}_bold_space-MNI152NLin2009cAsym_variant-smoothAROMAnonaggr_preproc.nii.gz
    MASK=$WD/${SUBJECT}/func/${SUBJECT}_task-${TASK}rest_run-${RUN}_bold_space-MNI152NLin2009cAsym_brainmask.nii.gz
    OUTPUT=$OUTDIR/${TASK}/$SUBJECT/zFC_150-${RUN}.csv
    CONFOUNDS=$WD/${SUBJECT}/func/${SUBJECT}_task-${TASK}rest_run-${RUN}_bold_confounds.tsv
    GLOBAL=$WD/${SUBJECT}/func/${SUBJECT}_task-${TASK}rest_run-${RUN}_GlobalSignal.tsv
    FRAMEWISE=$WD/${SUBJECT}/func/${SUBJECT}_task-${TASK}rest_run-${RUN}_FD.tsv

    if [ -f "$INPUT" ]; 
    then
      # extract only relevant columns
      awk '{print $3}' $CONFOUNDS > $GLOBAL
      awk '{print $7}' $CONFOUNDS > $FRAMEWISE
      FD=`./extract_connectivity.py $INPUT $ATLAS $MASK $CONFOUNDS $OUTPUT 2.0 full 0`
#      ./extract_connectivity.py $INPUT $ATLAS $MASK $CONFOUNDS $OUTPUT 2.0 full 0
      echo "${SUBJECT};${TASK};${RUN};$FD" >> $OUTDIR/FramewiseDisplacement.csv

      # select based on criteria
      FD_MEAN=`echo $FD | cut -d';' -f1`
      FD_MAX=`echo $FD | cut -d';' -f2`
      FD_3=`echo $FD | cut -d';' -f4`
      INVALID=`echo $FD | cut -d';' -f5`

    if [ $(bc <<< "$FD_MEAN <= $MEAN_THR") -eq 1 ] && 
       [ $(bc <<< "$FD_MAX <= $MAX_THR") -eq 1 ] && 
       [ $(bc <<< "$FD_3 <= $PROP_THR") -eq 1 ] &&
       [ $(bc <<< "$INVALID <= $MAX_INVALID") -eq 1 ];
    then
        LIST="$LIST $OUTPUT"
    fi

    fi
    done
    rm $OUTDIR/${TASK}/$SUBJECT/zFC_150.csv
    echo "*******************************************************************************"
    echo "LIST: $LIST"
    echo "*******************************************************************************"

    if [ "$LIST" != "" ]; 
    then
      echo $SUBJECT >> $OUTDIR/Subjects_${TASK}.txt
      ./average_connectivity.py $OUTDIR/${TASK}/$SUBJECT/zFC_150.csv "$LIST"
    fi
  done
done
