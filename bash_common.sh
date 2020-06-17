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
	echo -e "${RED}$MSG${NOCOLOR}"
	exit 1
}


function success
{
  echo ""
  echo -e "${GREEN}SUCCESS!${NOCOLOR}"
  exit 0
}