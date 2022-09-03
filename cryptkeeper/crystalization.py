"""
This file contains functions that utilize the output of more than one dependencies 
to derive new insights of the sequence of interest

:author: croots
"""

from collections import namedtuple

def build_unified_feature_list(expressed_ORFs, tss_predictions, rit_rhotermpredict_predictions, rit_transterm_predictions):
    """Combinds output of all dependencies to create a unified feature list"""
    feature_data = namedtuple('feature_data', ['feature_type', 'feature_start', 'feature_end', 'strength'])
    feature_list = []
    for feature in expressed_ORFs:
        feature_type = None
        feature_start = None
        feature_end = None
        strength = None
        feature_list.append(feature_data(feature_type, feature_start, feature_end, strength))
    for feature in tss_predictions:
        feature_type = None
        feature_start = None
        feature_end = None
        strength = None
        feature_list.append(feature_data(feature_type, feature_start, feature_end, strength))
    for feature in rit_rhotermpredict_predictions:
        feature_type = None
        feature_start = None
        feature_end = None
        strength = None
        feature_list.append(feature_data(feature_type, feature_start, feature_end, strength))
    for feature in rit_transterm_predictions:
        feature_type = None
        feature_start = None
        feature_end = None
        strength = None
        feature_list.append(feature_data(feature_type, feature_start, feature_end, strength))

    return feature_list

def identify_dsrna(feature_list):
    """Identifies dsRNA by looking for a region between two promoters with no other features"""
    pass

def identify_mutables():
    """Identifies regions that are known to be unstable within regions predicted to be expressed"""
    pass

