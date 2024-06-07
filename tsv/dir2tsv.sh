#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 3 ]; then
  echo "Combine all .tsv-files of a directory into one .tsv-file"
  echo "#1: directory with tab-delimited files"
  echo "#2: file with common header line if it is missing in files | ''"
  echo "#3: column name for file names"
  exit 1
fi
D=$1
H=$2
C="$3"


F=`ls $D | head -1 || true`
if [ -z "$F" ]; then
  error "Directory $D is empty"
fi
if [ ! -s $D/$F ]; then
  error "$D/$F is empty"
fi


TMP=`mktemp`


N=1
if [ $H ]; then
  cp $H $TMP
else
  head -1 $D/$F > $TMP
  N=2
fi
head -1 $TMP | grep -s '^#' > /dev/null || printf '#'
head -1 $TMP | sed 's/^/'"$C"'\t/1'                       

$THIS/../trav $D "tail -n +$N %d/%f | sed 's/^/'%f'\t/1'" 


rm $TMP
