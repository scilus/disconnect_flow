#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import glob
import json
import os
import pandas as pd
import re

def _build_arg_parser():
    parser = argparse.ArgumentParser(
        description='Extract and combine disconets',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('in_folder',
                        help='Results folder.')
    parser.add_argument('out_file',
                        help='csv file')
    add_overwrite_arg(parser)

    return parser


def main():
    parser = _build_arg_parser()
    args = parser.parse_args()
    assert_inputs_exist(parser, args.in_folder)
    assert_outputs_exist(parser, args, args.out_file)

    in_folder = args.in_folder
    print('Input folder : %s' % in_folder)

    atlases = ['CorticoCortical', 'CorticoStriatal', 'CorticoThalamic']
    indices = [[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17],
               [1,2,3,4,5,6,7,8,9,10,12,13,14,16,17],
               [1,2,3,4,5,6,7,8,9,10,11,13,15,16]]
    sides = ['L', 'R', 'LR', 'RL','']

    subjects = next(os.walk(in_folder))[1]

    stats = pd.DataFrame(columns=['subjectID', 'trkID'])


    expression = '(.*)_([0-9]{6}).*'
    check_expression = re.compile(expression)

    for idx, curr_atlas in enumerate(atlases):
        curr_stats = pd.DataFrame(columns=['subjectID', 'trkID'])
        for indice in indices[idx]:
            for side in sides:
                if side:
                    currColumn = ('_').join([curr_atlas, str(indice), side])
                else:
                    currColumn = ('_').join([curr_atlas, str(indice)])

                if side:
                    curr_search = os.path.join(in_folder,
                                           '*', # Subject Folder
                                           '*' + curr_atlas + '*', # Atlas Folder
                                           '*' + curr_atlas + '*_'+ str(indice) + '_' + side + '.txt')
                else:
                    curr_search = os.path.join(in_folder,
                                               '*', # Subject Folder
                                           '*' + curr_atlas + '*', # Atlas Folder
                                           '*' + curr_atlas + '*_'+ str(indice) + '.txt')



                curr_data = glob.glob(curr_search)
                curr_data.sort()
                for curr_json in curr_data:

                    check_re = check_expression.match(curr_json)
                    id = os.path.basename(check_re.groups()[0])
                    trk = check_re.groups()[1]
                    with open(curr_json) as f:
                        d = json.load(f)
                        if side:
                            sc_af_bdo = d['Filter_0']['streamline_count_after_filtering']
                        else:
                            sc_af_bdo = d['Filter_1']['streamline_count_after_filtering']
                        sc_af = d['streamline_count_final_filtering']


                    stats = add_stats(stats, id, trk, sc_af_bdo, sc_af, currColumn)
                    curr_stats = add_stats(curr_stats, id, trk, sc_af_bdo, sc_af, currColumn)

        print(curr_atlas)
        curr_stats.to_csv(curr_atlas+'.csv')

    stats.to_csv(args.out_file)

def add_stats(stats, id, trk, sc_af_bdo, sc_af, currColumn):

    for nColumn in [currColumn, currColumn + '_total']:
        if not nColumn in stats.columns:
            stats[nColumn] = ""

    index = stats[(stats['subjectID']==id) & (stats['trkID']==trk)].index
    val = ('/').join([str(sc_af), str(sc_af_bdo)])
    val = str(sc_af)
    val_total = str(sc_af_bdo)
    if index.empty:
        newEntry = pd.DataFrame([[id, trk, val]], columns=['subjectID', 'trkID', currColumn])
        stats = stats.append(newEntry, ignore_index=True)

        index = stats[(stats['subjectID']==id) & (stats['trkID']==trk)].index
        index = index[0]
        stats.loc[index][currColumn+'_total'] = val_total

    else:
        index = index[0]
        stats.loc[index][currColumn] = val
        stats.loc[index][currColumn+'_total'] = val_total

    return stats


def assert_inputs_exist(parser, required, optional=None):
    """
    Assert that all inputs exist. If not, print parser's usage and exit.
    :param parser: argparse.ArgumentParser object
    :param required: string or list of paths
    :param optional: string or list of paths.
                     Each element will be ignored if None
    """
    def check(path):
        if not os.path.isdir(path):
            parser.error('Input file {} does not exist'.format(path))

    if isinstance(required, str):
        required = [required]

    if isinstance(optional, str):
        optional = [optional]

    for required_file in required:
        check(required_file)
    for optional_file in optional or []:
        if optional_file is not None:
            check(optional_file)


def assert_outputs_exist(parser, args, required, optional=None):
    """
    Assert that all outputs don't exist or that if they exist, -f was used.
    If not, print parser's usage and exit.
    :param parser: argparse.ArgumentParser object
    :param args: argparse namespace
    :param required: string or list of paths
    :param optional: string or list of paths.
                     Each element will be ignored if None
    """
    def check(path):
        if os.path.isfile(path) and not args.overwrite:
            parser.error('Output file {} exists. Use -f to force '
                         'overwriting'.format(path))

    if isinstance(required, str):
        required = [required]

    if isinstance(optional, str):
        optional = [optional]

    for required_file in required:
        check(required_file)
    for optional_file in optional or []:
        if optional_file is not None:
            check(optional_file)

def add_overwrite_arg(parser):
    parser.add_argument(
        '-f', dest='overwrite', action='store_true',
        help='Force overwriting of the output files.')

if __name__ == "__main__":
    main()
