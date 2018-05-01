#!/bin/csh -f

if ($# != 1) then
  echo "#1: run qresub (0/1)"
  exit 1
endif


set N = 0
while (1)
  sleep 10  # PAR
  set Q = `qstat | grep -v '^job-ID' | grep -v '^---' | grep -v '   d[tr]   ' | head -1 | wc -l`
  if ($Q[1] == 0)  break
  
  @ N = $N + 1
  if ($N > 200) then  # PAR
    set N = 0
    set L = `qstat | grep -v '^job-ID' | grep -v '^---' | grep -v '   d[tr]   ' | grep '  [rT]  ' | sed 's/^ *//1' | cut -f 1 -d ' '`
    echo "Re-submitting $#L grid jobs ..."
    while ($#L)
      if ($1) then
	      qresub $L[1] -h u
	      if ($? == 0)  qdel -f $L[1]
	    else
	      qdel -f $L[1]
	    endif
      shift L
    end
    qrls  -h u  -u $USER > /dev/null
  endif
end


sync
