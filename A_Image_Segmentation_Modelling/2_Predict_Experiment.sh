#!/usr/bin/env bash
set -euo pipefail

# Run the trained Mask R-CNN on cropped images and their quadrants.
# Usage: bash 2_Predict_Experiment.sh <work_dir> <checkpoint_prefix>

if [[ $# -ne 2 ]]; then
    echo "Usage: bash $0 <work_dir> <checkpoint_prefix>" >&2
    exit 2
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORK_DIR="$(cd "$1" && pwd)"
CHECKPOINT_PREFIX="$(python -c 'import os,sys;print(os.path.abspath(sys.argv[1]))' "$2")"
PREDICT_SCRIPT="$SCRIPT_DIR/modelling_scripts/predict_json.py"

if [[ ! -f "$PREDICT_SCRIPT" ]]; then
    echo "Missing prediction script: $PREDICT_SCRIPT" >&2
    exit 1
fi
if [[ ! -d "$WORK_DIR/Border" || ! -d "$WORK_DIR/Split" ]]; then
    echo "Missing Border/ or Split/ directory. Run cropping first." >&2
    exit 1
fi
if [[ ! -f "${CHECKPOINT_PREFIX}.index" ]]; then
    echo "Missing TensorFlow checkpoint index: ${CHECKPOINT_PREFIX}.index" >&2
    exit 1
fi
shopt -s nullglob
checkpoint_data=("${CHECKPOINT_PREFIX}.data-"*)
if [[ ${#checkpoint_data[@]} -eq 0 ]]; then
    echo "Missing TensorFlow checkpoint data files: ${CHECKPOINT_PREFIX}.data-*" >&2
    exit 1
fi

# All images must have four quadrants; refuse mixed outputs from an earlier run.
images=0
for border_session in "$WORK_DIR"/Border/*; do
    [[ -d "$border_session" ]] || continue
    session="$(basename "$border_session")"
    split_session="$WORK_DIR/Split/$session"
    [[ -d "$split_session" ]] || { echo "Missing Split/$session" >&2; exit 1; }
    for border in "$border_session"/*.jpg; do
        [[ -f "$border" ]] || continue
        stem="$(basename "$border" .jpg)"
        for part in 0 1 2 3; do
            [[ -s "$split_session/${stem}_part${part}.jpg" ]] || {
                echo "Missing Split/$session/${stem}_part${part}.jpg" >&2
                exit 1
            }
        done
        ((images+=1))
    done
done
[[ "$images" -gt 0 ]] || { echo "No cropped images found." >&2; exit 1; }

existing=("$WORK_DIR"/Border/*/*.json "$WORK_DIR"/Split/*/*.json)
if [[ ${#existing[@]} -gt 0 ]]; then
    echo "Existing prediction JSON files detected; use a fresh working directory." >&2
    exit 1
fi

cd "$WORK_DIR"
for category in Split Border; do
    for session_dir in "$WORK_DIR/$category"/*; do
        [[ -d "$session_dir" ]] || continue
        session="$(basename "$session_dir")"
        inputs=()
        for file in "$session_dir"/*.jpg; do
            [[ -f "$file" ]] || continue
            inputs+=("$category/$session/$(basename "$file")")
        done
        [[ ${#inputs[@]} -gt 0 ]] || continue
        python "$PREDICT_SCRIPT" \
            --config 'MODE_FPN=True' \
            --output . \
            --load "$CHECKPOINT_PREFIX" \
            --predict "${inputs[@]}"
        for input in "${inputs[@]}"; do
            prediction="${input%.*}.json"
            if [[ ! -s "$WORK_DIR/$prediction" ]]; then
                echo "Prediction output missing: $WORK_DIR/$prediction" >&2
                exit 1
            fi
        done
    done
done

printf 'Prediction JSON files written beside cropped images and quadrants in %s\n' "$WORK_DIR"
