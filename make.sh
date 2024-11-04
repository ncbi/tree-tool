#!/bin/bash --noprofile
# For github distribution

unset AT_NCBI

cd "$( dirname $0 )"
export CPP_DIR=$PWD
source ./bash_common.sh
if [ $# -ne 1 ]; then
  echo "Make executables"
  echo "#1: rebuild (0/1)"
  exit 1
fi
REB=$1


function work
{
  local DIR=$1
  section $DIR
  cd $DIR
  if [ $REB == 1 ]; then
    rm *.o
  fi
  if [ -e Makefile ]; then
    make all
  fi
}


work $CPP_DIR
work $CPP_DIR/dm
work $CPP_DIR/dm/conversion
work $CPP_DIR/phylogeny
work $CPP_DIR/phylogeny/database
work $CPP_DIR/genetics
work $CPP_DIR/dissim/nw
work $CPP_DIR/dissim
work $CPP_DIR/ncbitax
work $CPP_DIR/tsv
work $CPP_DIR/web
work $CPP_DIR/xml


section "Software version"
$CPP_DIR/phylogeny/makeDistTree -version


success
