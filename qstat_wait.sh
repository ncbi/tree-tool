#!/bin/csh -f

while (1)
  sleep 10  # PAR
  set Q = `qstat | grep -v '^job-ID' | grep -v '^---' | grep -v '   d[tr]   ' | head -1 | wc -l`
  if ($Q[1] == 0)  break
end
