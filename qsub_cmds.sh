#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/bash_common.sh
N=100  # PAR
if [ $# -ne 4 ]; then
  echo "Submit a list of commands to the grid in the batches of $N and create #3.tsv"
  echo "#1: list of commands producing one-line output"
  echo "#2: intermediary input directory (not in /tmp, does not exist)"
  echo "#3: output directory (not in /tmp)"
  echo "#4: header file for the .tsv-output"
  exit 1
fi
CMDS=$1
IN=$2
OUT=$3
H=$3


TMP=$( mktemp )


mkdir $IN
$THIS/splitList $CMDS $N $IN -qc

mkdir $OUT
ls $IN > $TMP.in
while [ -s $TMP.in ] 
do
  $THIS/grid_wait.sh 1
  $THIS/trav $TMP.in "$QSUB_5 -N j%f %Qsource $IN/%f > $OUT/%f%Q > /dev/null"
  $THIS/qstat_wait.sh 360000 1  # PAR
  $THIS/trav $TMP.in "[ %D( < out/%f wc -l ) -eq %D( < cmd/%f wc -l ) ] || echo %f" > $TMP.in-new
  mv $TMP.in-new $TMP.in
done

rm -r $IN/ &

$THIS/tsv/dir2tsv.sh $OUT $H "" > $OUT.tsv
rm -r $OUT/ &


rm $TMP*
