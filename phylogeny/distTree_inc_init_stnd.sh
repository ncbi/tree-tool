#!/bin/bash
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 6 ]; then
  echo "Initialize an incremental distance tree directory with standard parameters for a specific biological project in $THIS/inc"
  echo "#1: output directory"
  echo "#2: biological project"
  echo "#3: SQL server name"
  echo "#4: database name on the SQL server"
  echo "#5: complete path to a directory for bulk insert into the database"
  echo "#6: path in Universal Naming Convention to the bulk directory #5"
  exit 1
fi
INC=$1
FROM=$2
SERVER=$3
DATABASE=$4
BULK_LOCAL=$5
BULK_REMOTE=$6


if [ ! -e $FROM/variance ]; then
  error "$FROM must be an incremental distance tree directory"
fi

( sqsh-ms  -S $SERVER  -D $DATABASE  -m bcp  -C "select @@version" > /dev/null ) || error "Database is not available"


$THIS/distTree_inc_init.sh $INC 1 "" 0 NAN NAN "" 0 $SERVER $DATABASE $BULK_LOCAL $BULK_REMOTE
cp $FROM/* $INC/


TMP=`mktemp`

CPP_DIR=$THIS/..
ls $FROM | grep '.sh$' > $TMP
$THIS/../trav $TMP -noprogress "sed %qs|CPP_DIR|"$CPP_DIR"|g%q $FROM/%f > $INC/%f"

rm $TMP*



