#!/usr/bin/env bash
# Tests, then the MIC outlier analysis (MHK and MGIT) -> output/mic/{mhk,mgit}/ (figures start from empty directories).
#SBATCH --job-name=kansasii-mic
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=/shares/sander.imm.uzh/MM/kansasii/output/mic/logs/%x-%j.log

set -euo pipefail
source "${KANSASII_ROOT:-/shares/sander.imm.uzh/MM/kansasii}/repos/kansasii-docs/templates/kansasii_paths.sh"
REPO=$REPOS/kansasii-mic
OUTPUT=$OUTPUT_ROOT/mic
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate env_mic
cd "$REPO"
export MPLBACKEND=Agg
pytest -q tests
rm -rf "$OUTPUT/mhk" "$OUTPUT/mgit"
python scripts/01_run_analysis.py --out-dir "$OUTPUT"

# provenance (see kansasii-docs/docs/PROVENANCE.md)
bash "$MANIFEST_SCRIPT" "$OUTPUT" "$REPO" submit_mic.sbatch "$MIC_CSV" "$SCREENING_MAP_LINK" "$SCREENING_MAP_AST_JC"
