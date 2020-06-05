#!/bin/bash
THIS=`dirname $0`
source $THIS/bash_common.sh
if [ $# -ne 2 ]; then
  echo "#1: seconds to wait after all jobs are in 'r' state"
  echo "#2: run qresub (0/1)"
  exit 1
fi
SECONDS=$1
QRESUB=$2


SLEEP_SEC=10  # PAR
PERIODS=$(( $SECONDS / $SLEEP_SEC ))


while true; do
  sleep $SLEEP_SEC
  $THIS/grid_wait.sh 0
  set +o errexit
  Q=`qstat | grep -v '^job-ID' | grep -v '^---' | grep '   qw   ' | head -1 | wc -l`
  set -o errexit
  if [ $Q == 0 ]; then
    break
  fi
done


N=0
while true; do
  sleep $SLEEP_SEC
  $THIS/grid_wait.sh 0
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
    $THIS/grid_wait.sh 1
    L=(`qstat | grep -v '^job-ID' | grep -v '^---' | grep -v '   d[tr]   ' | grep '  [rTt]  ' | sed 's/^ *//1' | cut -f 1 -d ' '`)
    set -o errexit
    if [ $QRESUB == 1 ]; then
      echo "Re-submitting ${#L[@]} grid jobs ..."
    else
      echo "Deleting ${#L[@]} grid jobs ..."
    fi
    i=0
    while [ $i -lt ${#L[@]} ]; do
      set +o errexit
      if [ $QRESUB == 1 ]; then
        $THIS/grid_wait.sh 1
	      qresub ${L[i]} -h u
	      S=$?
	      if [ $S == 0 ]; then
          $THIS/grid_wait.sh 1
	        qdel -f ${L[i]}
	      fi
	    else
        $THIS/grid_wait.sh 1
	      qdel -f ${L[i]}
	    fi
      set -o errexit
      i=$(( $i + 1 ))
    done
    if [ $QRESUB == 1 ]; then
      $THIS/grid_wait.sh 1
      qrls  -h u  -u $USER > /dev/null
    fi
  fi
done


sync
