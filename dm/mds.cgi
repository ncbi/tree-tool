#!/bin/bash --noprofile
source ./unCgi.sh


       NAME=`cat ${TMP}_name`
      CLASS=`cat ${TMP}_class`
DIST_LEGEND=`cat ${TMP}_distLegend`

DIST_ATTR=""
DIST_TYPE=""
if [ -e ${TMP}_distType ]; then
  DIST_ATTR=`cat ${TMP}_distAttr`
  DIST_TYPE=`cat ${TMP}_distType`
fi

mv ${TMP}_dm $TMP.dm


./head.sh "$NAME" white 0


if [ "$DIST_TYPE" == "" ]; then
  $PROD/cpp/dm/pca                        -maxClusters 6  -class "$CLASS"  -json $TMP.json  -verbose 0  -mds                      $TMP pca >& $TMP.mds
else 
  $PROD/cpp/dm/mds  -attrType $DIST_TYPE  -maxClusters 6  -class "$CLASS"  -json $TMP.json  -verbose 0        -attr "$DIST_ATTR"  $TMP     >& $TMP.mds
fi

echo "<script>"
echo "var mds = "
  cat $TMP.json
  echo ";"
  echo "var distLegend = '$DIST_LEGEND';"
  echo 'var pageName = "'$NAME'";'
  cat ./pcplot.js
echo "</script>"


source ./foot.sh
