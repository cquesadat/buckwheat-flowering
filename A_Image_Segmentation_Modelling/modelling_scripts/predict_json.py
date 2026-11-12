#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import itertools
import numpy as np
import os
import tensorflow as tf
import cv2
import tqdm
import json



import tensorpack.utils.viz as tpviz
from tensorpack.predict import MultiTowerOfflinePredictor, OfflinePredictor, PredictConfig
from tensorpack.tfutils import SmartInit, get_tf_version_tuple
from tensorpack.tfutils.export import ModelExporter
from tensorpack.utils import fs, logger
from skimage.measure import find_contours

from source.data_utils.dataset import DatasetRegistry
from source.data_utils.buckwheat import register_buckwheat
from source.utils.config import config as cfg
from source.utils.config import finalize_configs
from source.data_utils.data import get_eval_dataflow, get_train_dataflow
from eval import DetectionResult, multithread_predict_dataflow, predict_image
from source.modelling.generalized_rcnn import ResNetC4Model, ResNetFPNModel
from source.utils.viz import (
    draw_annotation, draw_final_outputs, draw_predictions,
    draw_proposal_recall, draw_final_outputs_blackwhite)


def do_visualize(model, model_path, nr_visualize=100, output_dir='output'):
    """
    Visualize some intermediate results (proposals, raw predictions) inside the pipeline.
    """
    df = get_train_dataflow()
    df.reset_state()

    pred = OfflinePredictor(PredictConfig(
        model=model,
        session_init=SmartInit(model_path),
        input_names=['image', 'gt_boxes', 'gt_labels'],
        output_names=[
            'generate_{}_proposals/boxes'.format('fpn' if cfg.MODE_FPN else 'rpn'),
            'generate_{}_proposals/scores'.format('fpn' if cfg.MODE_FPN else 'rpn'),
            'fastrcnn_all_scores',
            'output/boxes',
            'output/scores',
            'output/labels',
        ]))

    if os.path.exists(output_dir):
        raise FileExistsError('Visualization output directory already exists: {}'.format(output_dir))
    fs.mkdir_p(output_dir)
    with tqdm.tqdm(total=nr_visualize) as pbar:
        for idx, dp in itertools.islice(enumerate(df), nr_visualize):
            img, gt_boxes, gt_labels = dp['image'], dp['gt_boxes'], dp['gt_labels']

            rpn_boxes, rpn_scores, all_scores, \
                final_boxes, final_scores, final_labels = pred(img, gt_boxes, gt_labels)

            # draw groundtruth boxes
            gt_viz = draw_annotation(img, gt_boxes, gt_labels)
            # draw best proposals for each groundtruth, to show recall
            proposal_viz, good_proposals_ind = draw_proposal_recall(img, rpn_boxes, rpn_scores, gt_boxes)
            # draw the scores for the above proposals
            score_viz = draw_predictions(img, rpn_boxes[good_proposals_ind], all_scores[good_proposals_ind])

            results = [DetectionResult(*args) for args in
                       zip(final_boxes, final_scores, final_labels,
                           [None] * len(final_labels))]
            final_viz = draw_final_outputs(img, results)

            viz = tpviz.stack_patches([
                gt_viz, proposal_viz,
                score_viz, final_viz], 2, 2)

            if os.environ.get('DISPLAY', None):
                tpviz.interactive_imshow(viz)
            cv2.imwrite("{}/{:03d}.png".format(output_dir, idx), viz)
            pbar.update()


def do_evaluate(pred_config, output_file):
    num_tower = max(cfg.TRAIN.NUM_GPUS, 1)
    graph_funcs = MultiTowerOfflinePredictor(
        pred_config, list(range(num_tower))).get_predictors()

    for dataset in cfg.DATA.VAL:
        logger.info("Evaluating {} ...".format(dataset))
        dataflows = [
            get_eval_dataflow(dataset, shard=k, num_shards=num_tower)
            for k in range(num_tower)]
        all_results = multithread_predict_dataflow(dataflows, graph_funcs)
        output = output_file + '-' + dataset
        DatasetRegistry.get(dataset).eval_inference_results(all_results, output)

def save_co_ordinates(image, results):
    image = os.path.basename(image)
    image_data = []

    for r in results:
        mask = r.mask
        padded_mask = np.zeros(
            (mask.shape[0] + 2, mask.shape[1] + 2), dtype=np.uint8)
        padded_mask[1:-1, 1:-1] = mask
        contours = find_contours(padded_mask, 0.5)
        for verts in contours:
            verts = np.fliplr(verts) - 1
            list_co_ordinates = np.moveaxis(verts, 1, 0).tolist()

            region = {"shape_attributes": {"name": "polygon", 
                                           "all_points_x": list_co_ordinates[0],
                                           "all_points_y": list_co_ordinates[1]},
                      "region_attributes": {}}
            image_data.append(region)
    #size = img.shape[0]*img.shape[1]
    data = {"filename": image, "regions": image_data}
    #json_file = {image:data}
    return data

  


def do_predict(pred_func, input_file):
    img = cv2.imread(input_file, cv2.IMREAD_COLOR)
    if img is None:
        raise ValueError("Unable to read image: {}".format(input_file))

    results = predict_image(img, pred_func)
    coord_data = save_co_ordinates(input_file, results)
    image_size = img.shape[0] * img.shape[1]
    coord_data["size"] = image_size
    json_file = {coord_data["filename"] + str(image_size): coord_data}

    # Preserve relative input directories below the selected output directory.
    if os.path.isabs(input_file):
        json_name = os.path.splitext(os.path.basename(input_file))[0] + ".json"
    else:
        json_name = os.path.splitext(os.path.normpath(input_file))[0] + ".json"
        if json_name == '..' or json_name.startswith('..' + os.sep):
            raise ValueError('Input paths must not escape the working directory')
    json_output = os.path.join(args.output, json_name)
    os.makedirs(os.path.dirname(os.path.abspath(json_output)), exist_ok=True)
    with open(json_output, "w") as outfile:
        json.dump(json_file, outfile)

    logger.info("Segmentation JSON written to {}".format(json_output))



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--load', help='load a model for evaluation.', required=True)
    parser.add_argument('--visualize', action='store_true', help='visualize intermediate results')
    parser.add_argument('--evaluate', help="Run evaluation. "
                                           "This argument is the path to the output json evaluation file")
    parser.add_argument('--predict', help="Run prediction on a given image. "
                                          "This argument is the path to the input image file", nargs='+')

    parser.add_argument('--output', help="Output path of the prediction. ", default='output/')
    parser.add_argument('--benchmark', action='store_true', help="Benchmark the speed of the model + postprocessing")
    parser.add_argument('--config', help="A list of KEY=VALUE to overwrite those defined in config.py",
                        nargs='+')
    parser.add_argument('--output-pb', help='Save a model to .pb')
    parser.add_argument('--output-serving', help='Save a model to serving file')

    args = parser.parse_args()
    if args.config:
        cfg.update_args(args.config)

    register_buckwheat(cfg.DATA.BASEDIR)

    MODEL = ResNetFPNModel() if cfg.MODE_FPN else ResNetC4Model()

    if not tf.test.is_gpu_available():
        from tensorflow.python.framework import test_util
        assert get_tf_version_tuple() >= (1, 7) and test_util.IsMklEnabled(), \
            "Inference requires either GPU support or MKL support!"
    assert args.load
    finalize_configs(is_training=False)

    if args.predict or args.visualize:
        cfg.TEST.RESULT_SCORE_THRESH = cfg.TEST.RESULT_SCORE_THRESH_VIS

    if args.visualize:
        do_visualize(MODEL, args.load)
    else:
        predcfg = PredictConfig(
            model=MODEL,
            session_init=SmartInit(args.load),
            input_names=MODEL.get_inference_tensor_names()[0],
            output_names=MODEL.get_inference_tensor_names()[1])

        if args.output_pb:
            ModelExporter(predcfg).export_compact(args.output_pb, optimize=False)
        elif args.output_serving:
            ModelExporter(predcfg).export_serving(args.output_serving)

        if args.predict:
            predictor = OfflinePredictor(predcfg)
            for image_file in args.predict:
                do_predict(predictor, image_file)
        elif args.evaluate:
            assert args.evaluate.endswith('.json'), args.evaluate
            do_evaluate(predcfg, args.evaluate)
        elif args.benchmark:
            df = get_eval_dataflow(cfg.DATA.VAL[0])
            df.reset_state()
            predictor = OfflinePredictor(predcfg)
            for _, img in enumerate(tqdm.tqdm(df, total=len(df), smoothing=0.5)):
                # This includes post-processing time, which is done on CPU and not optimized
                # To exclude it, modify `predict_image`.
                predict_image(img[0], predictor)
