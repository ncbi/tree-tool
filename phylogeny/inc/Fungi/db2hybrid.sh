#!/bin/csh -f

if ($# != 1) then
  echo "Input: uniColl..Hybrid"
  echo "#1: Genome.id"
  exit 1
endif

sqsh-ms  -S PROTEUS  -U anyone  -D uniColl  -P allowed  -L exit_failcount=1 << EOF 
  EXEC Genome_hybrid $1, 0.1  -- PAR
  go -m bcp
EOF
if ($?) exit 1
