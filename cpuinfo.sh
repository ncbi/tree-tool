#!/bin/bash
THIS=`dirname $0`
source $THIS/bash_common.sh

grep MHz /proc/cpuinfo | sed 's/^.*: //1' | $THIS/dm/count
