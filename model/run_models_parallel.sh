#!/bin/bash
# run several models with different initializations if needed

INPUT_FILE=$1
SAMPLES_DIR=$2
STAN_FILE=$3
PRIORS_FILE=$4
FIT_FC=$5
USEMCMC=$6
TOTAL_INITS=$7
TOTAL_CHAINS=4

mkdir $SAMPLES_DIR
rm $SAMPLES_DIR/vbsamples*reduced $SAMPLES_DIR/vbsamples*simplified $SAMPLES_DIR/vbsamples #$SAMPLES_DIR/vbsamples_?.csv
doit(){
if [ $USEMCMC -eq 1 ]; then
Rscript fit_models_connectome.R $INPUT_FILE $SAMPLES_DIR/vbsamples $STAN_FILE $PRIORS_FILE $FIT_FC $USEMCMC > $SAMPLES_DIR/log
for FILE in $SAMPLES_DIR/vbsamples_?.csv; do
myname=`basename $FILE .csv | sed 's/_//'`
echo $myname
./reduce_samples.py $FILE $SAMPLES_DIR/${myname}_reduced muc sigmac raw -mcmc > $SAMPLES_DIR/log
./simplify_samples.py $FILE $SAMPLES_DIR/${myname}_simplified muc sigmac raw >> $SAMPLES_DIR/log
#rm $FILE 
done
rm $SAMPLES_DIR/log

else

Rscript fit_models_connectome.R $INPUT_FILE $SAMPLES_DIR/vbsamples${1} $STAN_FILE $PRIORS_FILE $FIT_FC $USEMCMC > $SAMPLES_DIR/log${1}
./reduce_samples.py $SAMPLES_DIR/vbsamples${1} $SAMPLES_DIR/vbsamples${1}_reduced muc sigmac raw > $SAMPLES_DIR/log${1}

CONVERGED=`cat $SAMPLES_DIR/log${1} | grep CONVERGED |wc -l`
if [ $CONVERGED==1 ]; then 
  rm $SAMPLES_DIR/log${1} #$SAMPLES_DIR/vbsamples${1}*; 
fi
fi
}

MAXJOBS=10
for i in `seq $TOTAL_INITS`; do

PROCESSES=`ps | grep run_models |wc -l`

while [ $PROCESSES -gt $MAXJOBS ]; do
       
       PROCESSES=`ps | grep run_models |wc -l`
       sleep 30                    
done
doit $i 

done

if [ $USEMCMC -eq 1 ]; then
while [ `echo $SAMPLES_DIR/vbsamples*_simplified | wc -w` -lt $TOTAL_CHAINS ]; do  sleep 300; done

# collect samples
cp $SAMPLES_DIR/vbsamples1_reduced $SAMPLES_DIR/vbsamples
for i in `seq 2 $TOTAL_CHAINS`; do
 if [ `cat $SAMPLES_DIR/vbsamples${i}_reduced |wc -l` -gt 1 ]; then 
	sed '1d' $SAMPLES_DIR/vbsamples${i}_reduced >> $SAMPLES_DIR/vbsamples
 fi

done
#rm $SAMPLES_DIR/vbsamples*_reduced


else

while [ `echo $SAMPLES_DIR/vbsamples*_reduced | wc -w` -lt $TOTAL_INITS ]; do  sleep 300; done

# collect samples in one file
cp $SAMPLES_DIR/vbsamples1_reduced $SAMPLES_DIR/vbsamples
for i in `seq 2 $TOTAL_INITS`; do
 if [ `cat $SAMPLES_DIR/vbsamples${i}_reduced |wc -l` -gt 1 ]; then 
	sed '1d' $SAMPLES_DIR/vbsamples${i}_reduced >> $SAMPLES_DIR/vbsamples
 fi

done

rm $SAMPLES_DIR/vbsamples*_reduced

fi

