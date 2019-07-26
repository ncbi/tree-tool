#/bin/bash
source bash_common.sh

sed -i 's|^CPP_DIR =.*$|CPP_DIR = '$PWD'|1' MakeRules

make all

cd dm
make all

cd ../phylogeny
make all
