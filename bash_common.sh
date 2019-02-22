set -o nounset
set -o errexit
set -o posix
set -o pipefail
export LC_ALL=C

QSUB_="qsub  -r n  -cwd  -V  -P unified  -j y  -v SGE_NOMAIL  -v SGE_SUMMARY=/dev/null  -b y   -o /dev/null  -l h_vmem=36G,mem_free=500M,m_mem_free=500M,reserve_mem=500M"
export QSUB_5="$QSUB_,h_rt=36000" 
