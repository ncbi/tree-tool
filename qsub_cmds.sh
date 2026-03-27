#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/bash_common.sh
N=100  # PAR
if [ $# -ne 3 ]; then
  echo "Submit a list of commands to the grid in the batches of $N"
  echo "#1: list of commands"
  echo "#2: intermediary input directory (not in /tmp, does not exist)"
  echo "#3: output directory (not in /tmp)"
  exit 1
fi
CMDS=$1
IN=$2
OUT=$3


mkdir $IN
$THIS/splitList $CMDS $N $IN -qc

mkdir $OUT
$THIS/grid_wait.sh 1
$THIS/trav $IN "$QSUB_5 -N j%f %Qsource %d/%f > $OUT/%f%Q > /dev/null"
$THIS/qstat_wait.sh 360000 1  # PAR


