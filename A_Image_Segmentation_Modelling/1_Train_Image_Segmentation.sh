#!/usr/bin/env bash
set -euo pipefail

# Train the segmentation model using distinct annotated train and val splits.
# Usage: bash 1_Train_Image_Segmentation.sh <training_data_dir> <coco_weights.npz> <log_dir> [train.py options]

if [[ $# -lt 3 ]]; then
    echo "Usage: bash $0 <training_data_dir> <coco_weights.npz> <log_dir> [train.py options]" >&2
    exit 2
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_BASEDIR="$(cd "$1" && pwd)"
COCO_WEIGHTS="$2"
LOGDIR="$3"
shift 3

for split in train val; do
    if [[ ! -f "$DATA_BASEDIR/$split/merged.json" ]]; then
        echo "Missing $DATA_BASEDIR/$split/merged.json" >&2
        exit 1
    fi
done
if [[ ! -f "$COCO_WEIGHTS" ]]; then
    echo "Missing COCO initialization weights: $COCO_WEIGHTS" >&2
    exit 1
fi
if [[ -e "$LOGDIR" ]]; then
    echo "Training log directory already exists. Supply a new path: $LOGDIR" >&2
    exit 1
fi

python "$SCRIPT_DIR/modelling_scripts/train.py" \
    --config "DATA.BASEDIR=$DATA_BASEDIR" \
    "DATA.VAL=('val',)" "DATA.TRAIN=('train',)" \
    'MODE_FPN=True' 'TRAIN.LR_SCHEDULE=(1200,)' \
    --logdir "$LOGDIR" --load "$COCO_WEIGHTS" "$@"
