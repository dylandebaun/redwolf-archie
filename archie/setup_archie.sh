#!/bin/bash
# setup_archie.sh — fetch ArchIE at the pinned upstream commit, apply the
# local modifications used in this project (see patches/), and build ms.
#
#   patches/calc_stats_window_data.patch  real data: S* positions stay aligned
#                                         with the sites (after the reference
#                                         filter AND after dropping sites where
#                                         the focal hap is ancestral)
#   patches/calc_stats_ms_rw.patch        simulations: adds calc_stats_ms_rw.py
#                                         (red-wolf-major labels, same S* fixes)
#   msmodified/ms                         compiled from the unmodified upstream
#                                         source (the committed binary is macOS)
#
#   bash archie/setup_archie.sh            # -> archie/ArchIE/
#
# Safe to re-run: ArchIE is reset to the pinned upstream commit before the
# patches are applied, so updated patches always replace older ones.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
COMMIT=00820bc61beca8a958f04fb17291ff66e03b9e94   # sriramlab/ArchIE HEAD, 2026-10
if [[ ! -d ArchIE ]]; then
  git clone https://github.com/sriramlab/ArchIE.git ArchIE
fi
git -C ArchIE checkout --quiet "${COMMIT}"
git -C ArchIE checkout --quiet -- .     # drop previously applied patches
git -C ArchIE clean -q -f -d            # ... and files they added
for p in patches/*.patch; do
  git -C ArchIE apply "../${p}"         # fails loudly if a patch doesn't apply
  echo "applied ${p}"
done

# build the modified ms (same command as upstream create_training.sh)
# -std=gnu89: ms is K&R C; newer gcc/clang defaults reject implicit declarations
gcc -std=gnu89 -Wno-implicit-function-declaration -Wno-return-type -o ArchIE/msmodified/ms ArchIE/msmodified/ms.c ArchIE/msmodified/streec.c ArchIE/msmodified/rand1.c -lm -w
echo "built ArchIE/msmodified/ms"
