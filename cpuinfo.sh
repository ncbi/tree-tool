#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/bash_common.sh

grep -w "MHz" /proc/cpuinfo | sed 's/^.*: //1' | $THIS/count
