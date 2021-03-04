#!/bin/bash --noprofile
# For github distribution
cd `dirname $0`
source ./bash_common.sh


unset AT_NCBI
export CPP_DIR=$PWD

make all

cd dm
make all

cd ../phylogeny
make all

cd ../genetics
make all

cd ../dissim
make all

cd ..


echo ""
echo "Software version:"
phylogeny/makeDistTree -version


success
