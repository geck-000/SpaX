#!/bin/bash
# Verify a re-solve actually solved, rather than reporting success for nothing.
#
# Three ways this has silently failed already, each of which looks like a clean
# COMPLETED from Slurm:
#   - the skip guard finds a truncated ODB and exits 0 in seconds
#   - the solve is killed and leaves an ODB that opens as corrupt
#   - the extraction globs a load case that was never solved and writes an
#     empty CSV
# So check elapsed time, ODB freshness, and populated rows -- not job state.
set -u
W=/scratch/project_2019020/test_rve
R="ssh roihu"

echo "=== solve states ==="
$R "sacct -j 671167,671169 --format=JobID%14,State,Elapsed -X -n | grep -v '^$'" 2>/dev/null

echo
echo "=== jobs that exited suspiciously fast (< 60s) ==="
$R "sacct -j 671167,671169 --format=JobID%14,Elapsed -X -n -P | awk -F'|' '{split(\$2,t,\":\"); s=t[1]*3600+t[2]*60+t[3]; if(s<60) print \$1, \$2}'" 2>/dev/null

echo
echo "=== skip-guard hits ==="
$R "cd $W && grep -h 'SKIP:' logs/*671167* logs/*671169* 2>/dev/null | wc -l" 2>/dev/null

echo
echo "=== ODBs written today ==="
$R "cd $W && echo \"  ERGL -ben fresh: \$(find . -maxdepth 1 -name 'Job-ERGL_*-ben.odb' -newermt today | wc -l) of \$(ls Job-ERGL_*-ben.odb 2>/dev/null | wc -l)\"; echo \"  TORL -tor fresh: \$(find . -maxdepth 1 -name 'Job-TORL_*-tor.odb' -newermt today | wc -l) of \$(ls Job-TORL_*-tor.odb 2>/dev/null | wc -l)\"" 2>/dev/null

echo
echo "=== populated results ==="
$R "cd $W && python3 spax_status.py" 2>/dev/null
