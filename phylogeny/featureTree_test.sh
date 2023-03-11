#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 1 ]; then
  echo "#1: go"
  echo "Time: 11 min."
  exit 1
fi


TMP=`mktemp`
comment $TMP
#set -x


DIR=$THIS/data/featureTree

mkdir $TMP.dir

cp $DIR/gene.tar.gz $TMP.dir/
cp $DIR/obj.list    $TMP.dir/

gunzip $TMP.dir/gene.tar.gz
tar  -xf $TMP.dir/gene.tar  -C $TMP.dir

$THIS/featureTree.sh $TMP.dir/obj $TMP.dir/gene

echo ""
$THIS/makeFeatureTree  -qc  -threads 10  -input_tree $TMP.dir/obj.tree  -features $TMP.dir/gene  -input_core $TMP.dir/obj.core  -use_time | grep -vw "^CHRON" | grep -v "^Tree from file:" > $TMP.dir/obj.featureTree
 
diff $TMP.dir/obj.core        $DIR/obj.core
diff $TMP.dir/obj.featureTree $DIR/obj.featureTree


rm -r $TMP*


success


