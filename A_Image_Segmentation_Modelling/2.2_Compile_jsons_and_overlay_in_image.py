#!/usr/bin/env python3
"""Merge Border/Split instance predictions and quantify flower cover per photograph."""

import argparse
import csv
import json
import math
import os
import re
import tempfile
from datetime import datetime
from pathlib import Path

from PIL import Image, ImageDraw
from shapely.geometry import Polygon
from shapely.ops import unary_union


IMAGE_NAME = re.compile(
    r'^(Plot[^_]+_[^_]+)_Session_(\d{4}-\d{2}-\d{2})\.jpe?g$', re.I
)
COLUMNS = ('Plot', 'ID', 'Date', 'Total_Area')


def read_design(path):
    with path.open(newline='', encoding='utf-8-sig') as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames or not {'Plot', 'FE_ID'}.issubset(reader.fieldnames):
            raise ValueError('Design must contain Plot and FE_ID columns: {}'.format(path))
        mapping = {}
        for row in reader:
            plot = row['Plot'].strip()
            genotype = row['FE_ID'].strip()
            if not plot or not genotype:
                raise ValueError('Empty plot or genotype in {}'.format(path))
            if plot in mapping:
                raise ValueError('Duplicate plot in design: {}'.format(plot))
            mapping[plot] = genotype
    return mapping


def key_from_image(name):
    match = IMAGE_NAME.fullmatch(name)
    if not match:
        raise ValueError('Unexpected image filename: {}'.format(name))
    plot, day = match.groups()[:2]
    datetime.strptime(day, '%Y-%m-%d')
    return plot, day


def read_prediction(path):
    if not path.is_file():
        raise FileNotFoundError('Missing predicted segmentation: {}'.format(path))
    with path.open(encoding='utf-8') as handle:
        data = json.load(handle)
    if isinstance(data, dict) and 'regions' in data:
        regions = data['regions']
    elif isinstance(data, dict) and len(data) == 1:
        entry = next(iter(data.values()))
        if not isinstance(entry, dict) or 'regions' not in entry:
            raise ValueError('Invalid prediction format: {}'.format(path))
        regions = entry['regions']
    else:
        raise ValueError('Invalid prediction format: {}'.format(path))
    if not isinstance(regions, list):
        raise ValueError('Prediction regions are not a list: {}'.format(path))
    return regions


def region_polygon(region, source):
    shape = region.get('shape_attributes', {})
    xs, ys = shape.get('all_points_x'), shape.get('all_points_y')
    if not isinstance(xs, list) or not isinstance(ys, list) or len(xs) != len(ys) or len(xs) < 3:
        raise ValueError('Invalid polygon coordinates: {}'.format(source))
    poly = Polygon(zip(xs, ys))
    if not poly.is_valid:
        raise ValueError('Self-intersecting or invalid polygon: {}'.format(source))
    return poly


def merged_predictions(border_path, split_paths, parts):
    polygons = []
    for region in read_prediction(border_path):
        polygons.append(region_polygon(region, border_path))

    with Image.open(parts[0]) as image:
        left_width, top_height = image.size
    with Image.open(parts[1]) as image:
        right_width, top_right_height = image.size
    with Image.open(parts[2]) as image:
        bottom_left_width, bottom_height = image.size
    with Image.open(parts[3]) as image:
        bottom_right_width, bottom_right_height = image.size
    if (top_height != top_right_height or left_width != bottom_left_width or
            right_width != bottom_right_width or bottom_height != bottom_right_height):
        raise ValueError('Inconsistent quadrant dimensions: {}'.format(parts[0]))

    offsets = ((0, 0), (left_width, 0), (0, top_height), (left_width, top_height))
    for predicted, (dx, dy) in zip(split_paths, offsets):
        for region in read_prediction(predicted):
            poly = region_polygon(region, predicted)
            # Shapely 1.x/2.x compatible translation without extra dependencies.
            translated = Polygon([(x + dx, y + dy) for x, y in poly.exterior.coords])
            polygons.append(translated)

    if not polygons:
        return [], 0.0

    merged = unary_union(polygons)
    if isinstance(merged, Polygon):
        shapes = [merged]
    elif merged.geom_type == 'MultiPolygon':
        shapes = list(merged.geoms)
    else:
        raise ValueError('Unexpected merged geometry: {}'.format(merged.geom_type))

    # Region polygons are the saved VIA-style segmentation representation.
    result = []
    for poly in shapes:
        coords = list(poly.exterior.coords)
        result.append({
            'shape_attributes': {
                'name': 'polygon',
                'all_points_x': [x for x, _ in coords],
                'all_points_y': [y for _, y in coords]
            },
            'region_attributes': {}
        })

    # Area follows the exterior polygons represented in the exported JSONs.
    return result, sum(abs(Polygon(poly.exterior).area) for poly in shapes)


def write_output_json(path, filename, pixel_area, regions):
    content = {
        'filename': filename,
        'size': pixel_area,
        'regions': regions
    }
    with path.open('w', encoding='utf-8') as handle:
        json.dump({filename + str(pixel_area): content}, handle)


def write_overlay(image_path, output_path, regions, area):
    with Image.open(image_path) as source:
        image = source.convert('RGBA')
    overlay = Image.new('RGBA', image.size, (0, 0, 0, 0))
    draw = ImageDraw.Draw(overlay)
    for region in regions:
        attrs = region['shape_attributes']
        xy = list(zip(attrs['all_points_x'], attrs['all_points_y']))
        draw.polygon(xy, fill=(255, 0, 0, 80), outline=(255, 0, 0, 180))
    combined = Image.alpha_composite(image, overlay).convert('RGB')
    ImageDraw.Draw(combined).text((6, 6), 'Flower cover: {:.2f}%'.format(area),
                                  fill=(255, 255, 255))
    combined.save(output_path)


def inventory_originals(source_dir):
    original = source_dir / 'Original'
    if not original.is_dir():
        raise FileNotFoundError('Missing original photographs directory: {}'.format(original))
    expected = set()
    for session_dir in sorted(original.iterdir()):
        if not session_dir.is_dir():
            continue
        match = re.fullmatch(r'Session_(\d{2})_(\d{2})_(\d{4})', session_dir.name)
        if not match:
            raise ValueError('Unexpected session directory: {}'.format(session_dir))
        folder_day = '{}-{}-{}'.format(match.group(3), match.group(2), match.group(1))
        for photo in session_dir.iterdir():
            if not photo.is_file() or photo.suffix.lower() not in ('.jpg', '.jpeg'):
                continue
            key = key_from_image(photo.name)
            if key[1] != folder_day:
                raise ValueError('Image date differs from session: {}'.format(photo))
            if key in expected:
                raise ValueError('Duplicate source plot/date: {}'.format(key))
            expected.add(key)
    if not expected:
        raise ValueError('No original JPEG photographs found')
    return expected


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('experiment_dir', type=Path, help='Source experiment with Original/ and Design_<experiment>.csv')
    parser.add_argument('work_dir', type=Path, help='Working directory with Border/ and Split/')
    parser.add_argument('--overlays', action='store_true', help='Also generate annotated PNGs')
    args = parser.parse_args()

    source = args.experiment_dir.resolve()
    work = args.work_dir.resolve()
    name = source.name
    design_path = source / 'Design_{}.csv'.format(name)
    design = read_design(design_path)
    expected = inventory_originals(source)

    border_root = work / 'Border'
    split_root = work / 'Split'
    if not border_root.is_dir() or not split_root.is_dir():
        raise FileNotFoundError('Border/ and Split/ are required in {}'.format(work))
    output_csv = work / 'Results_Flowering_Perc_{}.csv'.format(name)
    if output_csv.exists():
        raise FileExistsError('Output CSV already exists: {}'.format(output_csv))

    # Generated files remain inside the working directory.
    records = []
    seen = set()
    pred_root = work / 'Predicted_images'
    for border_session in sorted(border_root.iterdir()):
        if not border_session.is_dir():
            continue
        split_session = split_root / border_session.name
        if not split_session.is_dir():
            raise FileNotFoundError('Missing split session: {}'.format(split_session))
        for border_img in sorted(border_session.glob('*.jpg')):
            key = key_from_image(border_img.name)
            plot, day = key
            if key not in expected:
                raise ValueError('Cropped photograph not in source experiment: {}'.format(border_img))
            if key in seen:
                raise ValueError('Duplicate plot/date: {}'.format(key))
            if plot not in design:
                raise ValueError('Plot not found in experiment design: {}'.format(plot))

            stem = border_img.stem
            split_images = [split_session / '{}_part{}.jpg'.format(stem, n) for n in range(4)]
            split_jsons = [p.with_suffix('.json') for p in split_images]
            for path in split_images:
                if not path.is_file():
                    raise FileNotFoundError('Missing quadrant: {}'.format(path))
            border_json = border_img.with_suffix('.json')
            regions, area_px = merged_predictions(border_json, split_jsons, split_images)

            with Image.open(border_img) as image:
                width, height = image.size
            if not width or not height:
                raise ValueError('Invalid cropped image dimensions: {}'.format(border_img))
            total_area = 100.0 * area_px / (width * height)
            if not math.isfinite(total_area) or total_area < 0 or total_area > 100:
                raise ValueError('Invalid flower-cover percentage for {}: {}'.format(border_img, total_area))
            dest = pred_root / border_img.parent.name
            dest.mkdir(parents=True, exist_ok=True)
            write_output_json(dest / (stem + '.json'), border_img.name, width * height, regions)
            if args.overlays:
                write_overlay(border_img, dest / (stem + '.png'), regions, total_area)
            records.append((plot, design[plot], day, total_area))
            seen.add(key)

    if seen != expected:
        missing = sorted(expected - seen)
        raise ValueError('{} source photographs lack complete predictions; examples: {}'.format(len(missing), missing[:10]))
    if not records:
        raise ValueError('No complete predictions were found')

    # Write once, after all rows have passed the completeness checks.
    with tempfile.NamedTemporaryFile(mode='w', newline='', encoding='utf-8',
                                     prefix='.flower_cover_', suffix='.csv',
                                     dir=str(work), delete=False) as handle:
        tmp = Path(handle.name)
        writer = csv.writer(handle, lineterminator='\n')
        writer.writerow(COLUMNS)
        for plot, genotype, day, area in records:
            writer.writerow((plot, genotype, day, area))
    try:
        os.replace(str(tmp), str(output_csv))
    finally:
        if tmp.exists():
            tmp.unlink()

    print('Processed {} photographs'.format(len(records)))
    print('Flower-cover CSV: {}'.format(output_csv))
    print('Merged predictions: {}'.format(pred_root))


if __name__ == '__main__':
    main()
