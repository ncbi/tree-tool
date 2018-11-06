#!/bin/bash
source bash_common.sh
if [ $# -ne 2 ]; then
  echo "Restore #1/search/#2/"
  echo "#1: Directory containing search/"
  echo "#2: New object"
  exit 1
fi


cp /dev/null $1/new/$2
rm -r $1/search/$2

