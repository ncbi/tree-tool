#!/bin/bash
THIS=`dirname $0`
source $THIS/bash_common.sh
if [ $# -ne 1 ]; then
  echo "Wait until grid is available"
  echo "#1: go"
  exit 1
fi


N=0
while [ 1 == 1 ]; do
  set +o errexit
  qstat &> /dev/null
  S=$?
  set -o errexit
  if [ $S == 0 ]; then
    break
  fi
  sleep 60
  N=$(( $N + 1 ))
done

if [ $N -gt 0 ]; then
  echo "GRID was waited for $N min."
fi
