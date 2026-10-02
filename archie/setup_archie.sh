#!/bin/bash
# setup_archie.sh — fetch ArchIE at the pinned upstream commit, apply the
# local modifications used in this project (see patches/), and build ms.
#
#   patches/calc_stats_window_data.patch  real data: S* uses the positions of
#                                         the reference-filtered sites
#   patches/calc_stats_ms_rw.patch        simulations: adds calc_stats_ms_rw.py
#                                         (red-wolf-major labels, same S* fix)
#   msmodified/ms                         compiled from the unmodified upstream
#                                         source (the committed binary is macOS)
#
#   bash archie/setup_archie.sh            # -> archie/ArchIE/
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
COMMIT=00820bc61beca8a958f04fb17291ff66e03b9e94   # sriramlab/ArchIE HEAD, 2026-10
if [[ ! -d ArchIE ]]; then
  git clone https://github.com/sriramlab/ArchIE.git ArchIE
fi
git -C ArchIE checkout --quiet "${COMMIT}"
for p in patches/*.patch; do
  if git -C ArchIE apply --check "../${p}" 2>/dev/null; then
    git -C ArchIE apply "../${p}" && echo "applied ${p}"
  else
    echo "skip ${p} (already applied or does not apply cleanly)"
  fi
done

# build the modified ms (same command as upstream create_training.sh)
gcc -o ArchIE/msmodified/ms ArchIE/msmodified/ms.c ArchIE/msmodified/streec.c ArchIE/msmodified/rand1.c -lm -w
echo "built ArchIE/msmodified/ms"
