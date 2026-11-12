#!/usr/bin/env bash
set -euo pipefail

# Prepare cropped images and four quadrants in a new working directory.
# Usage: bash 0_Crop_Images_Experiment.sh <experiment_dir> <work_dir>

if [[ $# -ne 2 ]]; then
    echo "Usage: bash $0 <experiment_dir> <work_dir>" >&2
    exit 2
fi

EXPERIMENT_DIR="$(cd "$1" && pwd)"
ORIGINAL_DIR="${EXPERIMENT_DIR}/Original"
WORK_DIR="$2"

if [[ ! -d "$ORIGINAL_DIR" ]]; then
    echo "Missing original photographs directory: $ORIGINAL_DIR" >&2
    exit 1
fi
if [[ -e "$WORK_DIR" ]]; then
    echo "Working directory already exists; use a new empty path: $WORK_DIR" >&2
    exit 1
fi

# Do not write generated files into the source experiment.
WORK_DIR="$(python -c 'import os,sys;print(os.path.abspath(sys.argv[1]))' "$WORK_DIR")"
case "${WORK_DIR}/" in
    "${EXPERIMENT_DIR}/"*)
        echo "Working directory must be outside the source experiment." >&2
        exit 1 ;;
esac

if command -v magick >/dev/null 2>&1; then
    IMAGE_TOOL=(magick)
elif command -v convert >/dev/null 2>&1; then
    IMAGE_TOOL=(convert)
else
    echo "ImageMagick is required (magick or convert)." >&2
    exit 1
fi

shopt -s nullglob nocaseglob nocasematch
sessions=("$ORIGINAL_DIR"/Session_*)
if [[ ${#sessions[@]} -eq 0 ]]; then
    echo "No Session_DD_MM_YYYY directories found under $ORIGINAL_DIR" >&2
    exit 1
fi

# Preflight every source image before creating output directories.
count=0
for session_dir in "${sessions[@]}"; do
    [[ -d "$session_dir" ]] || continue
    session="$(basename "$session_dir")"
    if [[ ! "$session" =~ ^Session_([0-9]{2})_([0-9]{2})_([0-9]{4})$ ]]; then
        echo "Unexpected session directory: $session" >&2
        exit 1
    fi
    session_date="${BASH_REMATCH[3]}-${BASH_REMATCH[2]}-${BASH_REMATCH[1]}"
    images=("$session_dir"/*.jpg "$session_dir"/*.jpeg)
    for image in "${images[@]}"; do
        [[ -f "$image" ]] || continue
        name="$(basename "$image")"
        if [[ ! "$name" =~ ^Plot[^_]+_[^_]+_Session_([0-9]{4}-[0-9]{2}-[0-9]{2})\.(jpg|jpeg)$ ]]; then
            echo "Unexpected photograph filename: $image" >&2
            exit 1
        fi
        if [[ "${BASH_REMATCH[1]}" != "$session_date" ]]; then
            echo "Photograph date does not match session: $image" >&2
            exit 1
        fi
        ((count+=1))
    done
done
if [[ "$count" -eq 0 ]]; then
    echo "No JPEG photographs found in the experiment." >&2
    exit 1
fi

mkdir -p "$WORK_DIR/Border" "$WORK_DIR/Split"
processed=0
for session_dir in "${sessions[@]}"; do
    [[ -d "$session_dir" ]] || continue
    session="$(basename "$session_dir")"
    images=("$session_dir"/*.jpg "$session_dir"/*.jpeg)
    [[ ${#images[@]} -gt 0 ]] || continue

    border_dir="$WORK_DIR/Border/$session"
    split_dir="$WORK_DIR/Split/$session"
    mkdir -p "$border_dir" "$split_dir"
    for image in "${images[@]}"; do
        [[ -f "$image" ]] || continue
        stem="$(basename "$image")"
        stem="${stem%.*}"
        cropped="$border_dir/$stem.jpg"
        "${IMAGE_TOOL[@]}" "$image" -gravity Center -crop 70%x70%+0+0 +repage "$cropped"
        "${IMAGE_TOOL[@]}" "$cropped" -crop 2x2@ +repage +adjoin "$split_dir/${stem}_part%d.jpg"
        for part in 0 1 2 3; do
            if [[ ! -s "$split_dir/${stem}_part${part}.jpg" ]]; then
                echo "Missing cropped quadrant for $image (part $part)" >&2
                exit 1
            fi
        done
        ((processed+=1))
    done
done

printf 'Prepared %s photographs (%s quadrant images) in %s\n' "$processed" "$((processed * 4))" "$WORK_DIR"
