#!/bin/bash
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 2 ]; then
  echo "Save the parameters and scripts of an incremental distance tree"
  echo "#1: Input incremental tree directory"
  echo "#2: Output directory"
  exit 1
fi
INC=$1
OUT=$2


rm -rf $OUT
mkdir $OUT
for F in `ls $INC/`; do
  if [ -d $INC/$F -o $F == "tree" -o $F == "dissim" ]; then
    continue;
  fi

  set +o errexit
  N=`file $INC/$F | grep -c 'ASCII text'`
  set -o errexit
  if [ $N -eq 0 ]; then
    continue
  fi
  
  N=`cat $INC/$F | wc -l`
  if [ $N -lt 20000 ]; then  # PAR
    set +o errexit
    EX=`file $INC/$F | grep -cw 'executable'`
    set -o errexit
    if [ $EX -eq 0 ]; then
      cp $INC/$F $OUT/$F
    else
      sed 's|/home/brovervv/code/cpp|CPP_DIR|g' $INC/$F | sed 's/-blastdb_version 4//g' > $OUT/$F
      chmod a+x $OUT/$F
    fi
  fi
done


rm -f $OUT/dissim
rm -f $OUT/dissim.bad
rm -f $OUT/outlier_genogroup
rm -f $OUT/tree
rm -f $OUT/version
rm -f $OUT/runlog

# May be confidential
rm -f $OUT/server
rm -f $OUT/database


