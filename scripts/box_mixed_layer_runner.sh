#!/bin/bash
# Study runner for the mixed-layer Krüger base-case campaign (vast.ai box).
# Queue format: TAG|EXTRA_FLAGS (one per line) in /root/study_queue.txt.
# Same retry/skip/wall-budget discipline as the K24-DEKNOB-1 runner.
set -u
cd /root/petch-current
export PETCH_DETERMINISTIC_EXCHANGE_WORKERS=12
COMMON="--dx-um 0.01 --radiosity-backend deterministic_extruded_2d --transport-device cuda:0 \
 --surface-state-remap-backend common_refinement --topology-change-policy continue_gas_cavity \
 --surface-model mixed_layer --max-wall-s 86400"
while IFS='|' read -r TAG EXTRA; do
  [ -z "${TAG:-}" ] && continue
  OUT=/root/petch-results/$TAG
  LOG=/root/petch-results/$TAG.log
  [ -f "$OUT/audit.json" ] && python -c "import json,sys;sys.exit(0 if json.load(open('$OUT/audit.json')).get('status')=='complete' else 1)" 2>/dev/null && { echo "SKIP_COMPLETE $TAG" >> /root/study.log; continue; }
  FAILS=0
  while true; do
    if [ -f "$OUT/checkpoint.npz" ]; then
      python scripts/krueger_2024_trench_pilot.py $COMMON $EXTRA --output $OUT --resume >> $LOG 2>&1
    else
      python scripts/krueger_2024_trench_pilot.py $COMMON $EXTRA --output $OUT >> $LOG 2>&1
    fi
    CODE=$?
    STATUS=$(python -c "import json;print(json.load(open('$OUT/audit.json')).get('status',''))" 2>/dev/null)
    if [ $CODE -eq 0 ] && [ "$STATUS" != "wall_budget_checkpoint" ]; then
      echo "RUN_COMPLETE $TAG status=$STATUS" >> /root/study.log; break
    fi
    if [ $CODE -eq 0 ]; then FAILS=0; continue; fi
    FAILS=$((FAILS+1)); echo "FAIL $TAG code=$CODE n=$FAILS" >> /root/study.log
    [ $FAILS -ge 3 ] && { echo "STOP $TAG" >> /root/study.log; break; }
    sleep 5
  done
done < /root/study_queue.txt
echo QUEUE_DONE >> /root/study.log
