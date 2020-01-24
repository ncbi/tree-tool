#!/bin/bash
source ~brovervv/code/cpp/bash_common.sh
if [ $# -ne 3 ]; then
  echo "Load a gzipped XML file"
  echo "#1: gzipped XML file"
  echo "#2: XML file number"
  echo "#3: schema"
  exit 1
fi
BIOS=$1
NUM=$2
SCH=$3


TMP=`mktemp`
echo $TMP


function load ()
#1: table = file
{
  echo $1
  bulk.sh $TMP.flat/$1 PathogenBiosample..$1 PATHOGEN_DETECT 0 0
}


gunzip $BIOS -c > $TMP
mkdir $TMP.flat
xml_schema2flat  $TMP $NUM $SCH $TMP.flat -qc

load BioSample                      
load BioSample_Attributes_Attribute 
load BioSample_Description_Synonym  
load BioSample_Ids_Id                   
load BioSample_Links_Link
load BioSample_Models_Model


rm -r $TMP*

