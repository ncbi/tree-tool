#! /bin/echo Please-source

QSUB_="qsub  -r n  -cwd  -V  -P unified  -j y  -v SGE_NOMAIL  -v SGE_SUMMARY=/dev/null  -b y   -o /dev/null"
export QSUB_5="$QSUB_  -l h_vmem=36G,mem_free=500M,m_mem_free=500M,h_rt=36000"   
export QSUB_L="$QSUB_  -l h_vmem=36G,mem_free=1G,m_mem_free=1G,h_rt=1000000"   

# osver=8 ??
