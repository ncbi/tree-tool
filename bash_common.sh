set -o nounset
set -o errexit
set -o posix
set -o pipefail
export LC_ALL=C

QSUB_="qsub  -r n  -cwd  -V  -P unified  -j y  -v SGE_NOMAIL  -v SGE_SUMMARY=/dev/null  -b y   -o /dev/null"
export QSUB_5="$QSUB_  -l h_vmem=36G,mem_free=500M,m_mem_free=500M,h_rt=36000"   # ,reserve_mem=500M
export QSUB_L="$QSUB_  -l h_vmem=36G,mem_free=1G,m_mem_free=1G,h_rt=360000"   # ,reserve_mem=1G

export NOCOLOR='\033[0m'
export RED='\033[1;31m'
export GREEN='\033[1;32m'
export YELLOW='\033[1;33m'


function error
{
	MSG="$1"
	echo -e "${RED}ERROR${NOCOLOR}"
	echo "$MSG"
	exit 1
}


function section
{
	MSG="$1"
  echo ""
  echo -e "${GREEN}${MSG}${NOCOLOR}"
}


function super_section
{
  MSG="$1"
  S=`echo '*' $MSG '*'`
  T=`echo "$S" | sed 's/./*/g'`
  echo ""
  section "$T\n$S\n$T"
}


function warning
{
	MSG="$1"
  echo ""
  echo -e "${YELLOW}${MSG}${NOCOLOR}"
}


function success
{
  section "SUCCESS!"
  exit 0
}
