#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/bash_common.sh
if [ $# -ne 6 ]; then
  echo "Submit a list of commands to the grid in the batches create #3.tsv"
  echo "#1: list of commands"
  echo "#2: intermediary input directory (not in /tmp, to be created)"
  echo "#3: output directory (not in /tmp, to be created)"
  echo "#4: batch size"
  echo "#5: header file for the .tsv-output"
  echo "#6: commands prodeuce one-line output (1/0)"
  exit 1
fi
CMDS=$1
IN=$2
OUT=$3
BATCH=$4
H=$5
ONE=$6


if [ -e $OUT.tsv ]; then
  error "$OUT.tsv exists"
fi


TMP=$( mktemp )
comment $TMP


mkdir $IN
sort -R $CMDS > $TMP.cmds
$THIS/splitList $TMP.cmds $BATCH $IN -qc

mkdir $OUT
ls $IN > $TMP.in
while [ -s $TMP.in ] 
do
  $THIS/grid_wait.sh 1
  $THIS/trav $TMP.in "$QSUB_L -N j%f %Qsource $IN/%f > $OUT/%f%Q > /dev/null"
  $THIS/qstat_wait.sh 360000 1  # PAR
  if [ $ONE -eq 0 ]; then
    break
  fi
  rm -f $TMP.in-new
  $THIS/trav $TMP.in "[ %D( < $IN/%f wc -l ) -eq %D( < $OUT/%f wc -l ) ]"  -errors $TMP.in-new
  mv $TMP.in-new $TMP.in
done

rm -r $IN/ &

$THIS/tsv/dir2tsv.sh $OUT $H "" > $OUT.tsv
rm -r $OUT/ &


rm $TMP*
wait


