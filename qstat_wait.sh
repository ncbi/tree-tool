#!/bin/csh -f

set N = 0
while (1)
  sleep 10  # PAR
  set Q = `qstat | grep -v '^job-ID' | grep -v '^---' | grep -v '   d[tr]   ' | head -1 | wc -l`
  if ($Q[1] == 0)  break
  
  @ N = $N + 1
  if ($N > 200) then  # PAR
    set N = 0
    set L = `qstat | grep -v '^job-ID' | grep -v '^---' | grep -v '   d[tr]   ' | grep '  r  ' | sed 's/^ *//1' | cut -f 1 -d ' '`
    echo "Re-submitting $#L grid jobs ..."
    if (! $#L)  exit 1
    while ($#L)
      qresub $L[1]
      qdel $L[1]
      shift L
    end
  endif
end

