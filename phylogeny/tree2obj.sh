#!/bin/bash
source bash_common.sh
if [ $# -ne 1 ]; then
  echo "Print objects of a distance tree, sorted"
  echo "#1: Distance tree"
  exit 1
fi

grep -v '^ *0x' $1 | sed 's/^ *//1' | sed 's/:.*$//1' | sort


