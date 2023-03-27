#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 2 ]; then
  echo "Combine all .tsv-files of a directory into one .tsv-file, save as #1.tsv"
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

if [ -e $D.tsv ]; then
  error "File $D.tsv exists"
fi

head -1 $D/$F | sed 's/$/\t'$C'/1'                       >  $D.tsv
$THIS/../trav $D "tail -n +2 %d/%f | sed 's/$/\t'%f'/1'" >> $D.tsv
