#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 3 ]; then
  echo "Check integrity of Tax.rank an dTax.sub_rank"
  echo "#1: Server name"
  echo "#2: Database name"
  echo "#3: Schema name"
  exit 1
fi
SERVER=$1
DB=$2
SCHEMA=$3


sqsh-ms  -S $SERVER  -D $DB << EOT 
  begin
    set nocount on

    print '$SCHEMA.Tax/Parent.rank:'
    select Tax.id, Tax.[rank], P.[rank] parent_rank
      from      $SCHEMA.Tax
           join $SCHEMA.Tax P on P.id = Tax.parent
      where Tax.[rank] < P.[rank]

    print ''
    print '$SCHEMA.Tax.sub_rank:'
    select C.id, C.[rank], C.sub_rank, P.sub_rank
      from      $SCHEMA.Tax C
           join $SCHEMA.Tax P on P.id = C.parent
      where     C.[rank] = P.[rank]
            and C.sub_rank != P.sub_rank + 1
  end
  go 
EOT





