#!/bin/bash --noprofile
THIS=$( dirname $0 )
source $THIS/../bash_common.sh
if [ $# -ne 2 ]; then
  echo "Save the parameters and scripts of an incremental distance tree"
  echo "#1: input incremental tree directory"
  echo "#2: output directory"
  exit 1
fi
INC=$1
OUT=$2


function bad_word
{
  FILE=$1
  WORD=$2  # 0/1
  WHAT=$3
  W=""
  if [ $WORD == 1 ]; then
    W="-w"
  fi
  if grep $W "$WHAT" $FILE; then
    error "File $FILE contains $WHAT"
  fi
}


rm -rf $OUT
mkdir -p $OUT
for F in $( ls $INC/ | grep -v "\.old$" ); do
  if [ -d $INC/$F -o $F == "tree" -o $F == "dissim" ]; then
    continue;
  fi

  if [ -s $INC/$F ]; then
    set +o errexit
    N=`file $INC/$F | grep -c 'ASCII text'`
    set -o errexit
    if [ $N -eq 0 ]; then
      continue
    fi
  fi
  
  N=`cat $INC/$F | wc -l`
  if [ $N -lt 20000 ]; then  # PAR
    set +o errexit
    EX=$( file $INC/$F | grep -cw 'executable' )
    set -o errexit
    if [ $EX -eq 0 ]; then
      cp $INC/$F $OUT/$F
    else
      bad_word $INC/$F 0 '~'
      bad_word $INC/$F 1 "home"
      bad_word $INC/$F 1 "panfs"
      bad_word $INC/$F 1 "netmnt"
      bad_word $INC/$F 1 "brovervv"
      bad_word $INC/$F 1 "anyone"
      bad_word $INC/$F 1 "allowed"
      bad_word $INC/$F 1 "PROTEUS"
      bad_word $INC/$F 1 "uniColl"
      bad_word $INC/$F 1 "LIST"
      bad_word $INC/$F 1 "LISTC"
      bad_word $INC/$F 1 "loadLIST"
      bad_word $INC/$F 1 "loadLISTC"
      bad_word $INC/$F 1 "marker"
      bad_word $INC/$F 1 "GenBank"
      L=$( grep '[^]}];' $INC/$F | grep -v 'false;' | grep -v 'true;' | grep -vw "if" || true )  # ';' in SQL
      if [ "$L" ]; then
        error "$L"
      fi
      
      set +o errexit
      grep -w "PANFS" $INC/$F
      S=$?
      set -o errexit
      if [ $S -eq 0 ]; then
        error "PANFS is found in $INC/$F"
      fi
      
      sed 's|$BROVER_CPP|CPP_DIR|g' $INC/$F | sed 's/-blastdb_version 4//g' > $OUT/$F
      chmod a+x $OUT/$F
    fi
  fi
done


rm -f $OUT/dissim
rm -f $OUT/dissim.bad
rm -f $OUT/finished
rm -f $OUT/good.expanded
rm -f $OUT/hybrid.new
rm -f $OUT/indiscern
rm -f $OUT/leaf
rm -f $OUT/outlier-genogroup
rm -f $OUT/runlog
rm -f $OUT/search.list
rm -f $OUT/tree
rm -f $OUT/version

# May be confidential
rm -f $OUT/server
rm -f $OUT/database
rm -f $OUT/bulk_remote


