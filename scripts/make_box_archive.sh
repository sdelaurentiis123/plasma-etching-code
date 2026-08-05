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
)

git archive --prefix=petch-current/ "$REF" "${EXCLUDES[@]}" | gzip > "$OUT"

# Verify rather than assume: the tarball itself is the artifact that ships.
LEAKED="$(tar tzf "$OUT" | grep -icE 'partner|resona|challenge' || true)"
if [ "$LEAKED" != "0" ]; then
  echo "REFUSING: $LEAKED sensitive path(s) present in $OUT" >&2
  tar tzf "$OUT" | grep -iE 'partner|resona|challenge' >&2
  rm -f "$OUT"
  exit 1
fi

echo "archive OK: $OUT ($(du -h "$OUT" | cut -f1), $(tar tzf "$OUT" | wc -l | tr -d ' ') paths, 0 sensitive)"
