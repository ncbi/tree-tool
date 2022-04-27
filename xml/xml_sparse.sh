#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 1 ]; then
  echo "Print XML file with shorted lines"
  echo "#1: XML file"
  exit 1
fi
F=$1


sed 's/>/>\n/g' $F
