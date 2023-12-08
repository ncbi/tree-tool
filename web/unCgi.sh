#! /bin/echo Please-source
# Output: $TMP, $TMP-files, $Q, $PROD


#if [ "$HTTP_CONNECTION" != " " ; then  ??
echo "Content-type: text/html"
echo ""
#fi

export PROD=/netmnt/vast01/gp/home/brovervv/code
source $PROD/cpp/bash_common.sh
export PATH=/netopt/genbank/subtool/bin:/export/home/sybase/utils/bin:/usr/local/hmmer/3.1/bin:$PATH
TMP=`mktemp`
Q='"'

$PROD/cpp/web/unCgi $TMP -verbose 0



