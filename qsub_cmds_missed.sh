#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/bash_common.sh
if [ $# -ne 5 ]; then
  echo "Run missing commands"
  echo "#1: qsub_cmds.sh input directory"
  echo "#2: qsub_cmds.sh output directory"
  echo "#3: output directory"
  echo "#4: Batch #"
  echo "#5: use grid (0/1)"
  exit 1
fi
IN_ORIG=$1
OUT_ORIG=$2
OUT=$3
B=$4
GRID=$5


mkdir $OUT/$B
N=$( < $OUT_ORIG/$B  wc -l )
N=$(( N + 1 ))

CMD="%f > $OUT/$B/%n"
if [ $GRID -eq 1 ]; then
  CMD="$QSUB_5 -N j$B.%n %Q$CMD%Q > /dev/null"
fi

trav -start $N $IN_ORIG/$B "$CMD" -step 1


