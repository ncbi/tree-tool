#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/bash_common.sh
if [ $# -ne 5 ]; then
  echo "Bulk-insert"
  echo "#1: SQL server name"
  echo "#2: local bulk directory"
  echo "#3: path in Universal Naming Convention to the bulk directory"
  echo "#4: input '|'-delimited file"
  echo "#5: output table on server #1"
  exit 1
fi
SERVER=$1
BULK_LOCAL=$2
BULK_REMOTE=$3
IN=$4
TABLE=$5


#set -x


TMP=$( mktemp )
#comment $TMP 
rm $TMP
TMP=$( basename $TMP )

cp $IN $BULK_LOCAL/$TMP
unix2dos -o $BULK_LOCAL/$TMP &> /dev/null
chmod a+r $BULK_LOCAL/$TMP

sqsh-ms -S $SERVER  << EOT
  begin try
    set nocount on;
    delete from $TABLE;
    bulk insert $TABLE from '$BULK_REMOTE\\$TMP'
      with
      (
        BATCHSIZE = 100000
      , CHECK_CONSTRAINTS
      , FIRE_TRIGGERS 
      , fieldterminator = '|'
      );
    set nocount off;
  end try
  begin catch
    print error_message()
    raiserror ('error', 11, 1)
  end catch
  go -f -h
EOT

rm $BULK_LOCAL/$TMP




