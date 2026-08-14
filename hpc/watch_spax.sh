#!/bin/bash
# Emit one line per campaign as it reaches a full complement of populated cells,
# then exit. Also emits on two failures that would otherwise look like waiting:
# a drained queue with campaigns still short, and scratch filling up.
#
# The quota line matters more than it looks. At the limit Abaqus is killed
# mid-write and leaves a truncated ODB, so the campaign appears to have solved
# and the extraction fails later with "database file is corrupt" -- on whatever
# was writing at the time, not necessarily on the campaign that filled the disk.
# One crossing cost two campaigns' compute, so it is watched, not remembered.
prev=" "
warned=""
while true; do
  s=$(ssh roihu "cd /scratch/project_2019020/test_rve && python3 spax_status.py" 2>/dev/null)
  if [ -z "$s" ]; then sleep 600; continue; fi

  for f in $(echo "$s" | grep complete | awk '{print $1}'); do
    case "$prev" in
      *" $f "*) ;;
      *) echo "COMPLETE $f"; prev="$prev$f " ;;
    esac
  done

  qline=$(echo "$s" | grep '^QUOTA')
  case "$qline" in
    *CRITICAL*)
      if [ "$warned" != "critical" ]; then
        echo "QUOTA CRITICAL: $qline -- solves will write truncated ODBs"
        warned="critical"
      fi ;;
    *HIGH*)
      if [ -z "$warned" ]; then
        echo "QUOTA HIGH: $qline -- clear extracted ODBs before the next campaign"
        warned="high"
      fi ;;
  esac

  q=$(echo "$s" | awk '/^QUEUE/{print $2}')
  if echo "$s" | grep -q ALLDONE; then
    echo "ALLDONE: all six campaigns populated"
    break
  fi
  if [ "$q" = "0" ]; then
    echo "STALLED: slurm queue empty but campaigns still incomplete"
    echo "$s" | grep waiting
    break
  fi
  sleep 600
done
