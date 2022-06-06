#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/bash_common.sh
if [ $# -ne 2 ]; then
  echo "#1: seconds to wait after all jobs are in 'r' state"
  echo "#2: run qresub (0/1)"
  exit 1
fi
SECS=$1  # SECONDS is a system variable
QRESUB=$2


SLEEP_SEC=10  # PAR
PERIODS=$(( $SECS / $SLEEP_SEC ))


while true; do
  sleep $SLEEP_SEC
  $THIS/grid_wait.sh 0
  set +o errexit
  Q=`qstat | grep -v '^job-ID' | grep -v '^---' | grep '   qw   ' | wc -l`
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
  Q=`qstat | grep -v '^job-ID' | grep -v '^---' | grep -v '   d[tr]   ' | wc -l`
  set -o errexit
  if [ $Q == 0 ]; then
    break
  fi
  
  N=$(( $N + 1 ))
  if [ $N -gt $PERIODS ]; then  
    N=0
    $THIS/grid_wait.sh 1
    set +o errexit
    L=(`qstat | grep -v '^job-ID' | grep -v '^---' | grep -v '   d[tr]   ' | grep '  [rTt]  ' | sed 's/^ *//1' | cut -f 1 -d ' '`)
    set -o errexit
    M=${#L[@]}
    if [ $M -ne 0 ]; then 
      date
      if [ $QRESUB == 1 ]; then
        echo "Re-submitting $M grid jobs ..."
        for i in $M; do
          $THIS/grid_wait.sh 1
  	      qresub ${L[$i]} -h u || true
        done
      else
        echo "Deleting $M grid jobs ..."
      fi
      $THIS/grid_wait.sh 1
      qdel -f ${L[@]} || true
      if [ $QRESUB == 1 ]; then
        $THIS/grid_wait.sh 1
        qrls  -h u  -u $USER > /dev/null
      fi
    fi
  fi
  
  if [ $SECS == 0 ]; then
    break
  fi
done


#sync
