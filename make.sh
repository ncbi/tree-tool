#!/bin/bash
source ./bash_common.sh

sed 's|^CPP_DIR =.*$|CPP_DIR = '$PWD'|1' MakeRules > MakeRules.tmp
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
