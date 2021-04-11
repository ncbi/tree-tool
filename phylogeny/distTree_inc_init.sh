#!/bin/bash --noprofile
THIS=`dirname $0`
source $THIS/../bash_common.sh
if [ $# -ne 13 ]; then
  echo "Initialize an incremental distance tree directory"
  echo "#1: output directory"
  echo "#2: min. number of dissimilarity requests to use GRID (> 0)"
  echo "#3: variance: lin|linExp|..."
  echo "#4: hybridness_min (> 1); 0 - no hybrids"
  echo "#5: dissim_boundary (> 0 or NAN)"
  echo "#6: genogroup_barrier (> 0 or NAN)"
  echo "#7: complete path to a 'phenotype' directory or ''"
  echo "#8: large (0/1): 1 - files in #1/new/ and #7/ are grouped into subdirectories named file2hash(<file name>)"
  echo "#9: object2closest.sql (0/1): 1 - object2closest.sh queries an SQL database and the number of concurrent connections is restricted to 30"
  echo "#10: SQL server name or ''"
  echo "#11: database name on the SQL server or ''"
  echo "#12: complete path to a directory for bulk insert into the database or ''"
  echo "#13: path in Universal Naming Convention to the bulk directory #11 or ''"
  exit 1
fi
INC=$1
GRID_MIN=$2
VARIANCE="$3"
HYBRIDNESS_MIN=$4
DISSIM_BOUNDARY=$5
GENOGROUP_BARRIER=$6
PHEN="$7"
LARGE=$8
REQUEST_CLOSEST_SQL=$9
SERVER="${10}"
DATABASE="${11}"
BULK_LOCAL="${12}"
BULK_REMOTE="${13}"


if [ $GRID_MIN -le 0 ]; then
  error "Bad GRID_MIN"
fi

if [ -n "$SERVER" -a -z "$DATABASE" ]; then
  error "Database is empty"
fi

if [ -z "$SERVER" -a -n "$DATABASE" ]; then
  error "Server is empty"
fi

if [ -n "$BULK_LOCAL" -a -z "$DATABASE" ]; then
  error "Bulk directory with an empty database"
fi

if [ -n "$BULK_LOCAL" -a -z "$BULK_REMOTE" ]; then
  error "Local bulk directory with an empty remote bulk directory"
fi

if [ -z "$BULK_LOCAL" -a -n "$BULK_REMOTE" ]; then
  error "Remote bulk directory with an empty local bulk directory"
fi


mkdir $INC

touch $INC/dissim
touch $INC/leaf
touch $INC/tree
touch $INC/indiscern

mkdir $INC/new
mkdir $INC/search

echo "1" > $INC/version

mkdir $INC/hist

echo $GRID_MIN > $INC/object2closest.grid
echo $GRID_MIN > $INC/pairs2dissim.grid
touch $INC/runlog

echo "$VARIANCE"          > $INC/variance
echo $DISSIM_BOUNDARY     > $INC/dissim_boundary
echo $GENOGROUP_BARRIER   > $INC/genogroup_barrier
echo $HYBRIDNESS_MIN      > $INC/hybridness_min
echo $SERVER              > $INC/server
echo $DATABASE            > $INC/database

if [ $LARGE == 1 ]; then
  touch $INC/large
fi

if [ $REQUEST_CLOSEST_SQL == 1 ]; then
  touch $INC/object2closest.sql
fi



function create_script
{
	NAME=$1.sh
	#
	cp $THIS/distTree_inc/$NAME $INC/
#echo -e "${YELLOW}Implement $INC/$NAME !${NOCOLOR}"
}
create_script objects_in_tree
if [ $SERVER ]; then
  create_script outlier2db
fi
create_script pairs2dissim
create_script object2closest
create_script qc
create_script qc_object
if false; then  # deprecated
  if [ $HYBRIDNESS_MIN != 0 ]; then
  	create_script db2unhybrid
  	create_script hybrid2db
  fi
fi
if [ "$DISSIM_BOUNDARY" != "NAN" ]; then
  create_script genogroup2db
fi


if [ "$BULK_LOCAL" ]; then
  ln -s $BULK_LOCAL $INC/bulk
  echo $BULK_REMOTE > $INC/bulk_remote
fi

if [ "$PHEN" ]; then
  ln -s $PHEN $INC/phen
fi

if [ $LARGE == 1 ]; then
  $THIS/../trav 1000 -zero -start 0 "mkdir $INC/new/%n" -threads 10
fi

