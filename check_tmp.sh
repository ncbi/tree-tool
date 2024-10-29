#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/bash_common.sh

AVAILABLE=$( df -k /tmp | tail -1 | awk '{print $4};' )
if [ $AVAILABLE -lt 1000000 ]; then  # PAR
  error "Space available on /tmp is $AVAILABLE Kb"
fi
