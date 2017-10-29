#!/bin/csh -f

if ($# != 1) then
  echo "#1: go"
  exit 1
endif


set DIR = data/featureTree

cp $DIR/gene.tar.gz .
if ($?) exit 1
cp $DIR/obj.list .
if ($?) exit 1

gunzip gene.tar.gz
if ($?) exit 1
tar -xf gene.tar
if ($?) exit 1

featureTree.sh obj gene
if ($?) exit 1

diff obj.core $DIR/obj.core
if ($?) exit 1
diff obj.tree $DIR/obj.tree
if ($?) exit 1

rm obj.list
rm obj.tree
rm obj.core
rm -r gene/
rm gene.tar
