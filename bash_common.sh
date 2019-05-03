set -o nounset
set -o errexit
set -o posix
set -o pipefail
export LC_ALL=C

QSUB_="qsub  -r n  -cwd  -V  -P unified  -j y  -v SGE_NOMAIL  -v SGE_SUMMARY=/dev/null  -b y   -o /dev/null"
export QSUB_5="$QSUB_  -l h_vmem=36G,mem_free=500M,m_mem_free=500M,reserve_mem=500M,h_rt=36000" 
export QSUB_L="$QSUB_  -l h_vmem=36G,mem_free=1G,m_mem_free=1G,reserve_mem=1G,h_rt=360000" 
