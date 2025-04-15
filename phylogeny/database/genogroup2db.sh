#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../../bash_common.sh
if [ $# -ne 2 ]; then
  echo "Populate GenomeTaxroot.outlier for genonogroup outliers"
  echo "#1: inc/"
  echo "#2: File genogroup_table"
  exit 1
fi
INC=$1
TAB=$2


SERVER=$( cat $INC/server )
DATABASE=$( cat $INC/database )
BULK_REMOTE=$( cat $INC/bulk_remote )
TAXROOT=$( cat $INC/../taxroot )

$THIS/../../bulk.sh $SERVER $INC/bulk $BULK_REMOTE $TAB $DATABASE..ListC

sqsh-ms  -S $SERVER  -D $DATABASE  -L exit_failcount=1  -L bcp_rowsep="" << EOF
  EXEC Genogroup2outliers $TAXROOT, 0
  go -m bcp
EOF
