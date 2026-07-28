#!/usr/bin/env bash
# Preflight smoke of the EXACT production Krüger mixed-layer configuration.
#
# Runs a few real steps locally on CPU with every production flag set, so any
# integration friction (extrusion-guard trips, transport/chemistry contract
# breaks, refusals) surfaces on the laptop in minutes instead of on a rented
# GPU at minute 30. Run before every box deploy.
#
# Usage:  scripts/preflight_krueger.sh [duration_s]   (default 0.05)
set -u

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DURATION="${1:-0.05}"
OUT="$(mktemp -d)"
LOG="$OUT/preflight.log"

cd "$ROOT" || { echo "PREFLIGHT FAIL: cannot cd to repo root"; exit 1; }

echo "preflight: mixed_layer + grazing reflection, dx=0.01, ${DURATION}s, CPU"

timeout 240 python scripts/krueger_2024_trench_pilot.py \
    --surface-model mixed_layer \
    --grazing-ion-reflection literature_v1 \
    --dx-um 0.01 \
    --duration-s "$DURATION" \
    --transport-device cpu \
    --surface-state-remap-backend common_refinement \
    --topology-change-policy continue_gas_cavity \
    --output "$OUT/run" > "$LOG" 2>&1
CODE=$?

if [ "$CODE" -ne 0 ]; then
    echo "PREFLIGHT FAIL (exit $CODE). First error:"
    grep -aE "Error|raise|Traceback|Refus|guard" "$LOG" | head -3
    echo "  full log: $LOG"
    exit 1
fi

STATUS="$(python -c "import json,sys; print(json.load(open('$OUT/run/audit.json')).get('status',''))" 2>/dev/null)"
STEPS="$(grep -ac '^step=' "$LOG")"

case "$STATUS" in
    complete|running|wall_budget_checkpoint)
        echo "PREFLIGHT PASS: status=$STATUS, $STEPS steps clean"
        rm -rf "$OUT"
        exit 0
        ;;
    *)
        echo "PREFLIGHT FAIL: unexpected audit status '$STATUS' after $STEPS steps"
        echo "  full log: $LOG"
        exit 1
        ;;
esac
