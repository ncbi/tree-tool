#!/bin/bash
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 7 ]; then
  echo "Initialize an incremental distance tree directory with standard parameters for a specific biological project in $THIS/inc"
  echo "#1: output directory"
  echo "#2: biological project"
  echo "#3: complete path to a 'phenotype' directory or ''"
  echo "#4: SQL server name"
  echo "#5: database name on the SQL server"
  echo "#6: complete path to a directory for bulk insert into the database"
  echo "#7: path in Universal Naming Convention to the bulk directory #6"
  exit 1
fi
INC=$1
FROM=$2
PHEN="$3"
SERVER=$4
DATABASE=$5
BULK_LOCAL=$6
BULK_REMOTE=$7


if [ ! -e $FROM/variance ]; then
  error "$FROM must be an incremental distance tree directory"
fi


$THIS/phylogeny/distTree_inc_init.sh $INC 1 "" 0 NAN NAN $PHEN 0 $SERVER $DATABASE $BULK_LOCAL $BULK_REMOTE
cp $FROM/* $INC/


