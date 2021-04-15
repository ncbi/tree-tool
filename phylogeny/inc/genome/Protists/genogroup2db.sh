#!/bin/bash --noprofile
source CPP_DIR/bash_common.sh
if [ $# -ne 1 ]; then
  echo "$0"
  echo "#1: File genogroup_table"
  exit 1
fi


loadLISTC $1

sqsh-ms  -S PROTEUS  -D uniColl  -L exit_failcount=1 << EOF | sed 's/|$//1' 
  EXEC Genogroup2outliers 2759 /*PAR*/, 0;
  go -m bcp
EOF

