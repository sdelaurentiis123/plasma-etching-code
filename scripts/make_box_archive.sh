#!/usr/bin/env bash
# Build the tarball shipped to rented compute, with the partner-sensitive
# exclusion codified rather than retyped.
#
# Why this exists: the exclusion was carried in operating instructions as the
# single pattern ':!PARTNER_ETCH_CHALLENGE*'.  A second sensitive file
# (RESONA_PATTERN_TRANSFER_PARTNERSHIP_*) is also tracked, matched no pattern,
# and therefore shipped to every box of the campaign until ml23 caught it.
# One pattern per file is a footgun; this script owns the list and verifies the
# result instead of trusting the invocation.
#
# 2026-08-06: a path-only control is not sufficient.  Three tracked roadmap
# docs NAME the partner in their body text and shipped to every box because
# their filenames are innocuous.  The control is now two-layer: paths are
# excluded, and the packed CONTENT is scanned so a future doc that mentions a
# partner cannot ship silently.  Content hits fail the build closed.
#
# Usage: scripts/make_box_archive.sh [OUT_TGZ] [GIT_REF]

set -euo pipefail

OUT="${1:-/tmp/petch-box.tgz}"
REF="${2:-HEAD}"

# Every path that must never reach rented hardware.  Add here, nowhere else.
EXCLUDES=(
  ':!PARTNER_*'
  ':!RESONA_*'
  ':!*PARTNERSHIP*'
  ':!*_CHALLENGE_*'
  # Body-text mentions found by the 2026-08-06 content scan.  None is needed
  # on a box; all are planning/orchestration docs.  The content scan below is
  # the backstop -- if a future doc names a partner, the build fails closed
  # rather than relying on this list being complete.
  ':!COMPETITIVE_DIFFERENTIATION_ROADMAP_*'
  ':!PHYSICS_AI_ACCELERATION_ROADMAP_*'
  ':!PHYSICS_FIRST_UNIFIED_ENGINE_*'
  ':!FRONTIER_LOOP.md'
  ':!ORCHESTRATION.md'
  ':!RESEARCH_DIFFERENTIABLE_TRANSPORT_*'
  ':!UNIFIED_ENGINE_VALIDATION_EXECUTION_PROGRAM_*'
  ':!VALIDATION_FIRST_SUPERSET_CAMPAIGN_*'
)

# Tokens that must not appear in shipped file CONTENT, not just in path names.
# Word-boundary matched: bare "Resona" is the partner, "resonance"/"resonant"/
# "resonator" are ordinary physics and must not trip the gate.
CONTENT_TOKENS='(^|[^a-z])resona([^a-z]|$)|partner_etch'


git archive --prefix=petch-current/ "$REF" "${EXCLUDES[@]}" | gzip > "$OUT"

# Verify rather than assume: the tarball itself is the artifact that ships.
LEAKED="$(tar tzf "$OUT" | grep -icE 'partner|resona|challenge' || true)"
if [ "$LEAKED" != "0" ]; then
  echo "REFUSING: $LEAKED sensitive path(s) present in $OUT" >&2
  tar tzf "$OUT" | grep -iE 'partner|resona|challenge' >&2
  rm -f "$OUT"
  exit 1
fi

# Second layer: scan what is actually inside the shipped files.
CONTENT_HITS="$(tar xzOf "$OUT" 2>/dev/null | grep -icE "$CONTENT_TOKENS" || true)"
# This script carries the tokens by construction; it is not a leak.
SELF_HITS="$(tar xzOf "$OUT" petch-current/scripts/make_box_archive.sh 2>/dev/null | grep -icE "$CONTENT_TOKENS" || true)"
CONTENT_HITS=$(( CONTENT_HITS - SELF_HITS ))
if [ "$CONTENT_HITS" != "0" ]; then
  echo "REFUSING: $CONTENT_HITS partner mention(s) in shipped file CONTENT" >&2
  echo "  offending files:" >&2
  for f in $(tar tzf "$OUT" | grep -E '\.(md|py|sh|txt|json)$' | grep -v make_box_archive); do
    if tar xzOf "$OUT" "$f" 2>/dev/null | grep -qiE "$CONTENT_TOKENS"; then
      echo "    $f" >&2
    fi
  done
  rm -f "$OUT"
  exit 1
fi

echo "archive OK: $OUT ($(du -h "$OUT" | cut -f1), $(tar tzf "$OUT" | wc -l | tr -d ' ') paths, 0 sensitive paths, 0 content mentions)"
