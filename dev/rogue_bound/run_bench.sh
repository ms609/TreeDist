#!/bin/bash
# Interleaved A/B/C benchmark; each process: 1 warm-up + N timed runs.
cd /tmp/claude-0/-home-claude/37427f5a-6383-50c4-a6e2-4664447d3c94/scratchpad/bound
out=bench_results.txt
: > $out
for round in 1 2; do
  for w in rand30 near80 near150 mid80 mid150 rand80 rand150; do
    reps=5
    if [ "$w" = "rand150" ]; then reps=3; fi
    for l in lib_base lib_arr lib_new2; do
      Rscript bench.R $l $w $reps >> $out 2>&1
    done
  done
  for l in lib_base lib_arr lib_new2; do
    Rscript bench.R $l rand30 3 2 >> $out 2>&1
  done
done
echo DONE >> $out
