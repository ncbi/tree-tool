#/bin/bash

source bash_common.sh

make all

cd dm
make all

cd ../phylogeny
make all
