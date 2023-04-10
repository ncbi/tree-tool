#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/bash_common.sh
if [ $# -ne 1 ]; then
  echo "Wait until UGE is available"
  echo "#1: print message (0/1)"
  exit 1
fi
MSG=$1


if [ $MSG == 1 ]; then
  comment "Waiting for UGE"
fi

N=0
while true; do
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
  echo "UGE was waited for for $N min."
fi


#sync
