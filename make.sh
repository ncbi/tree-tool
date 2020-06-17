#!/bin/bash
# For github distribution
source ./bash_common.sh
unset AT_NCBI

if false; then
DIR=$PWD
sed 's|^CPP_DIR =.*$|CPP_DIR = '$DIR'|1' MakeRules > MakeRules.tmp
mv MakeRules.tmp MakeRules
fi

make all

cd dm
make all

cd ../phylogeny
make all

cd ../genetics
make all

cd ../dissim
make all


if false; then
cd ../phylogeny
set +o errexit
ls inc/*.sh       1>  tmp.list  2> /dev/null  
ls inc/*/*.sh     1>> tmp.list  2> /dev/null
ls inc/*/*/*.sh   1>> tmp.list  2> /dev/null
ls inc/*/*/*/*.sh 1>> tmp.list  2> /dev/null
set -o errexit
../trav tmp.list -noprogress "sed 's|CPP_DIR|'$DIR'|g' %f > %f.tmp && mv %f.tmp %f && chmod a+x %f"
rm tmp.list
fi


echo ""
echo "Software version:"
./makeDistTree -version

echo ""
echo -e ${GREEN}SUCCESS!${NOCOLOR}
