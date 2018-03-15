#!/bin/csh -f

if ($# != 3) then
  echo "#1: hash directory"
  echo "#2: hash file"
  echo "#3: hash type (uniColl..GenomeHash.type)"
  exit 1
endif


cat $1/$2 | sed 's/^/'$2'#'$3'#/1' | tr '#' '\t'
if ($?) exit 1

