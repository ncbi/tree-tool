#!/bin/bash --noprofile
# For github distribution

unset AT_NCBI

cd `dirname $0`
source ./bash_common.sh
export CPP_DIR=$PWD


make all

cd dm
make all

cd conversion
make all

cd ../../phylogeny
make all

cd ../genetics
make all

cd ../dissim
make all

cd ../tsv
make all

cd ../xml
make all

cd ..


section "Software version"
phylogeny/makeDistTree -version


success
