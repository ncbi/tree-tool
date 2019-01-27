#!/bin/bash
EXEC=~brovervv/code/cpp
source $EXEC/bash_common.sh
if [ $# -ne 2 ]; then
  echo "#1: seconds to wait"
  echo "#2: run qresub (0/1)"
  exit 1
fi
SECONDS=$1
QRESUB=$2


SLEEP_SEC=10  # PAR
PERIODS=$(( $SECONDS / $SLEEP_SEC ))


while  [ 1 == 1 ]; do
  sleep $SLEEP_SEC
  set +o errexit
  Q=`qstat | grep -v '^job-ID' | grep -v '^---' | grep '   qw   ' | head -1 | wc -l`
  set -o errexit
  if [ $Q == 0 ]; then
    break
  fi
done


N=0
while [ 1 = 1 ]; do
  sleep $SLEEP_SEC
  set +o errexit
  Q=`qstat | grep -v '^job-ID' | grep -v '^---' | grep -v '   d[tr]   ' | head -1 | wc -l`
  set -o errexit
  if [ $Q == 0 ]; then
    break
  fi
  
  N=$(( $N + 1 ))
  if [ $N -gt $PERIODS ]; then  
    N=0
    set +o errexit
    L=(`qstat | grep -v '^job-ID' | grep -v '^---' | grep -v '   d[tr]   ' | grep '  [rTt]  ' | sed 's/^ *//1' | cut -f 1 -d ' '`)
    set -o errexit
    echo "Re-submitting $#L grid jobs ..."
    while [ $#L -gt 0 ]; do
      if [ $QRESUB == 1 ]; then
        set +o errexit
	      qresub ${L[0]} -h u
	      S=$?
        set -o errexit
	      if [ $S == 0 ]; then
	        qdel -f ${L[0]}
	      fi
	    else
	      qdel -f ${L[0]}
	    fi
      shift L
    done
    qrls  -h u  -u $USER > /dev/null
  fi
done


sync
