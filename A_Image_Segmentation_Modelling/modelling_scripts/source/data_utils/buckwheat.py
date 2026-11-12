"""Buckwheat flower instance-segmentation dataset loader."""

import json
import os
from functools import partial

import numpy as np
from source.data_utils.dataset import DatasetSplit, DatasetRegistry

__all__ = ['register_buckwheat']


class BuckwheatDemo(DatasetSplit):
    def __init__(self, base_dir, split):
        if split not in ('train', 'val'):
            raise ValueError('Unknown dataset split: {}'.format(split))
        self.imgdir = os.path.join(os.path.expanduser(base_dir), split)
        if not os.path.isfile(os.path.join(self.imgdir, 'merged.json')):
            raise FileNotFoundError('Missing annotation file: {}'.format(
                os.path.join(self.imgdir, 'merged.json')))

    def load_roidbs(self, training=True):
        with open(os.path.join(self.imgdir, 'merged.json')) as handle:
            annotations = json.load(handle)

        if isinstance(annotations, dict):
            entries = list(annotations.values())
        elif isinstance(annotations, list):
            entries = annotations
        else:
            raise ValueError('merged.json must contain annotation records')

        records = []
        seen = set()
        for entry in entries:
            if isinstance(entry, dict) and 'filename' not in entry and len(entry) == 1:
                entry = next(iter(entry.values()))
            if not isinstance(entry, dict) or 'filename' not in entry:
                raise ValueError('Invalid annotation record in {}'.format(self.imgdir))

            filename = os.path.join(self.imgdir, entry['filename'])
            if filename in seen:
                raise ValueError('Duplicate image annotation: {}'.format(filename))
            if not os.path.isfile(filename):
                raise FileNotFoundError('Missing annotated image: {}'.format(filename))
            seen.add(filename)
            record = {'file_name': filename}

            if training:
                regions = entry.get('regions', [])
                if isinstance(regions, dict):
                    regions = list(regions.values())
                boxes, segmentations = [], []
                for region in regions:
                    shape = region.get('shape_attributes', {})
                    px = shape.get('all_points_x', [])
                    py = shape.get('all_points_y', [])
                    if len(px) != len(py) or len(px) < 3:
                        raise ValueError('Invalid annotation polygon: {}'.format(filename))
                    polygon = np.stack((px, py), axis=1).astype(np.float32) + 0.5
                    minxy = polygon.min(axis=0)
                    maxxy = polygon.max(axis=0)
                    boxes.append([minxy[0], minxy[1], maxxy[0], maxxy[1]])
                    segmentations.append([polygon])

                count = len(boxes)
                record['boxes'] = np.asarray(boxes, dtype=np.float32).reshape((-1, 4))
                record['segmentation'] = segmentations
                record['class'] = np.ones((count,), dtype=np.int32)
                record['is_crowd'] = np.zeros((count,), dtype=np.int8)
            else:
                record['image_id'] = filename

            records.append(record)
        return records

    def training_roidbs(self):
        return self.load_roidbs(training=True)

    def inference_roidbs(self):
        return self.load_roidbs(training=False)

    def eval_inference_results(self, results, output=None):
        if output is not None:
            with open(output, 'w') as handle:
                json.dump(results, handle)
        return results


def register_buckwheat(basedir):
    for split in ('train', 'val'):
        DatasetRegistry.register(split, partial(BuckwheatDemo, basedir, split))
        DatasetRegistry.register_metadata(split, 'class_names', ['BG', 'Flower'])
