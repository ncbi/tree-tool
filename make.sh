#!/bin/bash
source ./bash_common.sh
unset AT_NCBI

DIR=$PWD
sed 's|^CPP_DIR =.*$|CPP_DIR = '$DIR'|1' MakeRules > MakeRules.tmp
mv MakeRules.tmp MakeRules

make all

cd dm
make all

cd ../phylogeny
make all

cd ../genetics
make all

cd ../dissim
make all


cd ../phylogeny
set +o errexit
ls inc/*.sh       1>  tmp.list  2> /dev/null  
ls inc/*/*.sh     1>> tmp.list  2> /dev/null
ls inc/*/*/*.sh   1>> tmp.list  2> /dev/null
ls inc/*/*/*/*.sh 1>> tmp.list  2> /dev/null
set -o errexit
../trav tmp.list -step 1 "sed 's|CPP_DIR|'$DIR'|g' %f > %f.tmp && mv %f.tmp %f && chmod a+x %f"
rm tmp.list
