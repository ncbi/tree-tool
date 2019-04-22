#!/bin/bash
if [ $# -ne 1 ]; then
  echo "Record genogroup outliers in the database"
  echo "#1: File genogroup_table"
  exit 1
fi
IN=$1


exit 1