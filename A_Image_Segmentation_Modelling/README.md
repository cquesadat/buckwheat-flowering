# A - Flower image segmentation

This pipeline uses a Mask R-CNN model (ResNet-50/FPN) to identify buckwheat flowers in field photographs and calculate flower-cover percentage per image.

## Inputs

For each experiment (`ETH23`, `ETH24`, `LfL23`, `LfL24`, `ETH24L`), provide:

- `Experiments/<EXPERIMENT>/Original/Session_DD_MM_YYYY/*.jpg` — field photographs named `Plot<ROW>_<RANGE>_Session_YYYY-MM-DD.jpg`
- `Experiments/<EXPERIMENT>/Design_<EXPERIMENT>.csv` — plot-to-genotype mapping (`Plot`, `FE_ID`)
- `modelling_scripts/model-400000.index` and `model-400000.data-00000-of-00001` — trained model checkpoint

The supplied `Results_Flowering_Perc_<EXPERIMENT>.csv` files are the archived flower-cover measurements; they can be used directly by Pipeline B.

## Run inference

Run from `A_Image_Segmentation_Modelling/`, setting a **new, nonexistent** working directory outside the experiment's input folder:

```bash
EXP_DIR="$PWD/Experiments/ETH23"
WORK_DIR="/path/to/new_work_directory/ETH23"
CHECKPOINT="$PWD/modelling_scripts/model-400000"

bash 0_Crop_Images_Experiment.sh "$EXP_DIR" "$WORK_DIR"
bash 2_Predict_Experiment.sh "$WORK_DIR" "$CHECKPOINT"
python 2.2_Compile_jsons_and_overlay_in_image.py "$EXP_DIR" "$WORK_DIR"
```

The crop step prepares the central 70% of each image and four quadrants. Prediction writes JSON polygons; the final step merges them and computes flower cover. Add `--overlays` to the final command to save annotated PNGs.

## Outputs

The working directory contains `Border/` (central crops), `Split/` (quadrants), `Predicted_images/` (merged predictions) and `Results_Flowering_Perc_<EXPERIMENT>.csv` with columns **`Plot,ID,Date,Total_Area`**. `Total_Area` is the percentage of cropped image area covered by flowers, with one row per original photograph.

## Optional model training

Provide `modelling_scripts/modelling_data/train/` and `val/`, each containing images and `merged.json`, plus `modelling_scripts/COCO-MaskRCNN-R50FPN2x.npz`. From this directory run:

```bash
bash 1_Train_Image_Segmentation.sh \
  "$PWD/modelling_scripts/modelling_data" \
  "$PWD/modelling_scripts/COCO-MaskRCNN-R50FPN2x.npz" \
  "/path/to/new_training_output"
```

Requires Bash, ImageMagick and a compatible Python 3.7 / TensorFlow 1.15.5 / TensorPack environment with the dependencies in `modelling_scripts/setup.py`. Upstream code attribution is retained in `modelling_scripts/source/utils/README.md`.
