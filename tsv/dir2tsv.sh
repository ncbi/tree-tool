#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 2 ]; then
  echo "Combine all .tsv-files of a directory into one .tsv-file"
  echo "#1: directory"
  echo "#2: column name for file names"
  exit 1
fi
D=$1
C="$2"


F=`ls $D | head -1 || true`
if [ -z "$F" ]; then
  error "Directory $D is empty"
fi

head -1 $D/$F | sed 's/$/\t'$C'/1'
$THIS/../trav $D "sed 's/$/\t'%f'/1' %d/%f" | grep -v '^#'
