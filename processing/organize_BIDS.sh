#!/bin/sh

#organize the data according to BIDS standard
ANATDIR='/home/benjamin.garzon/Data/DAD/data/T1Img'
FILES_DIR=/home/benjamin.garzon/Software/DAD/DAandConnectomes/processing

WD='/home/benjamin.garzon/Data/DAD/data_BIDS'
rm -r $WD
mkdir $WD
cd $ANATDIR

echo "{\"Name\": \"DAD\", \"BIDSVersion\": \"1.0.2\"}" > $WD/dataset_description.json
echo "participant_id" > $WD/participants.tsv


RUN=1

for NAME in D*; do
echo $NAME
echo $NAME >> $WD/participants.tsv

mkdir $WD/sub-${NAME}
mkdir $WD/sub-${NAME}/func/

  for TASK in GNG TAB RS; do

    json_text="{
    \"TaskName\": \"$TASK\",
    \"RepetitionTime\": 2.0,
    \"SliceTiming\": [0.0000,0.1081,0.2162,0.3243,0.4324,0.5405,0.6486,0.7568,0.8649,0.9730,1.0811,1.1892,1.2973,1.4054,1.5135,1.6216,1.7297,1.8378,1.9459,0.0541,0.1622,0.2703,0.3784,0.4865,0.5946,0.7027,0.8108,0.9189,1.0270,1.1351,1.2432,1.3514,1.4595,1.5676,1.6757,1.7838,1.8919],
    \"InstitutionName\": \"Karolinska Institute\"

    }"

    FMRIDIR="/home/benjamin.garzon/Data/DAD/data/FunImg_$TASK"  

    if [ "$TASK" == "RS" ]; then
    RUN=1
    cp $FMRIDIR/$NAME/data.nii $WD/sub-${NAME}/func/sub-${NAME}_task-${TASK}rest_run-${RUN}_bold.nii
    gzip $WD/sub-${NAME}/func/sub-${NAME}_task-${TASK}rest_run-${RUN}_bold.nii
    echo $json_text > $WD/sub-${NAME}/func/sub-${NAME}_task-${TASK}rest_run-${RUN}_bold.json
    fi

    if [ "$TASK" == "TAB" ]; then
    for RUN in 1 2; do
    cp $FMRIDIR/$NAME/TABrun${RUN}.nii.gz  $WD/sub-${NAME}/func/sub-${NAME}_task-${TASK}rest_run-${RUN}_bold.nii.gz
    echo $json_text > $WD/sub-${NAME}/func/sub-${NAME}_task-${TASK}rest_run-${RUN}_bold.json
    done
    fi

    if [ "$TASK" == "GNG" ]; then
    for RUN in 1 2 3; do
    cp $FMRIDIR/$NAME/GNGrun${RUN}.nii.gz $WD/sub-${NAME}/func/sub-${NAME}_task-${TASK}rest_run-${RUN}_bold.nii.gz
    echo $json_text > $WD/sub-${NAME}/func/sub-${NAME}_task-${TASK}rest_run-${RUN}_bold.json
    done
    fi

  done

  mkdir $WD/sub-${NAME}/anat/

  cp $ANATDIR/$NAME/coT1.nii $WD/sub-${NAME}/anat/sub-${NAME}_T1w.nii
  gzip $WD/sub-${NAME}/anat/sub-${NAME}_T1w.nii
done


