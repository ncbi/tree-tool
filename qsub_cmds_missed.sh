#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/bash_common.sh
if [ $# -ne 4 ]; then
  echo "Run missing commands"
  echo "#1: qsub_cmds.sh input directory"
  echo "#2: qsub_cmds.sh output directory"
  echo "#3: output directory"
  echo "#4: Batch #"
  exit 1
fi
IN_ORIG=$1
OUT_ORIG=$2
OUT=$3
B=$4


mkdir $OUT/$B
N=$( < $OUT_ORIG/$B  wc -l )
N=$(( N + 1 ))
trav -start $N $IN_ORIG/$B "$QSUB_5 -N j$B.%n %Q%f > $OUT/$B/%n%Q > /dev/null" -noprogress


