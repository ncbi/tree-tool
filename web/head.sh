#!/bin/bash --noprofile
if [ $# != 3 ]; then
  echo "#1: Title"
  echo "#2: bgcolor"
  echo "#3: DB is used (0/1)"
  echo "To be invoked by a cgi-script in the beginning after unCgi"
  exit 1
fi
TITLE=$1
BGCOLOR=$2
DB_USED=$3


echo '<html lang="en">'

# Internet Explorer
echo "<STYLE>"
echo "    hr { page-break-after: always }"
echo "</STYLE>"

echo "<head>"
  echo "<meta charset='utf-8' />"
  echo '<meta http-equiv="Cache-Control" content="no-store" />'
  echo "<title>$TITLE</title>"
  echo "<base href=https://intrawebdev8/staff/brovervv/>"
  if [ $DB_USED == 1 ]; then
    echo "<script type='text/javascript' src='jquery-1.11.1.js'></script>"
  fi
  echo '<script src="common.js?ver=1"></script>'
  if [ $DB_USED == 1 ]; then
    echo '<script src="common_db.js?ver=1"></script>'
  fi
echo "</head>"

echo "<body font='Lucida Console' bgcolor=$BGCOLOR onunload=''>"
echo "<h2><center>$TITLE</center></h2>"


