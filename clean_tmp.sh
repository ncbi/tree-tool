#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/bash_common.sh


TMP=$( mktemp )
comment $TMP


DIR=$( dirname $TMP )
NAME=$( basename $TMP )
ls -l $DIR | grep -w $USER | awk '{print $9};' | grep -vw $NAME > $TMP

trav $TMP "rm -r $DIR/%f"


rm $TMP*
