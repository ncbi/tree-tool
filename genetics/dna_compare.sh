#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../bash_common.sh
if [ $# -ne 2 ]; then
  echo "blastn 2 DNAs"
  echo "#1; DNA 1"
  echo "#2: DNA 2"
  exit 1
fi
F1=$1
F2=$2


TMP=$( mktemp )
#comment $TMP 
#set -x


cp $F2 $TMP
makeblastdb  -in $TMP   -dbtype nucl  -blastdb_version 4  -logfile /dev/null
blastn  -query $F1  -db $TMP  -dust no  


rm $TMP*
