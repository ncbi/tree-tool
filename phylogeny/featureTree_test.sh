#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 1 ]; then
  echo "#1: go"
  echo "Time: 11 min."
  exit 1
fi


DIR=data/featureTree

cp $DIR/gene.tar.gz .
cp $DIR/obj.list .

gunzip gene.tar.gz
tar -xf gene.tar

$THIS/featureTree.sh obj gene

echo ""
$THIS/makeFeatureTree  -qc  -threads 10  -input_tree obj.tree  -features gene  -input_core obj.core  -use_time | grep -vw "^CHRON" > obj.featureTree
 
diff obj.core $DIR/obj.core
diff obj.featureTree $DIR/obj.featureTree

rm obj.list
rm obj.tree
rm obj.core
rm obj.featureTree
rm -r gene/
rm gene.tar


success
