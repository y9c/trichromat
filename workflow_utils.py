#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2025 Ye Chang yech1990@gmail.com
# Distributed under terms of the GNU license.
#
# Created: 2025-05-10 03:18


import os
from collections import defaultdict


def resolve_path(path, workdir):
    if os.path.isabs(path):
        return path
    else:
        path = os.path.expanduser(path)
        if os.path.isabs(path):
            return path
        else:
            return os.path.relpath(path, workdir)


def _preprocess_reference(config, workdir):
    ref = defaultdict(list)
    for k, v in config["reference"].items():
        ref[k] = [resolve_path(v2, workdir) for v2 in v]
    if "genome" not in ref:
        ref["genome"] = []
    else:
        # check if genome records is more than 1, if so, raise error
        if len(ref["genome"]) > 1:
            raise ValueError(
                "Multiple genome records found in reference. Please provide only one genome record.\n"
                "Automatically join multiple genome records into one require rebuild the index file, which takes a lot of time.\n"
            )
    # if contamination in ref but is empty, remove it
    if "contamination" in ref and len(ref["contamination"]) == 0:
        del ref["contamination"]
    return ref


def _preprocess_reads(config, workdir):
    reads = defaultdict(lambda: defaultdict(dict))
    for s, v in config["samples"].items():
        for i, v2 in enumerate(v["data"], 1):
            reads[str(s)][f"run{i}"] = {
                k: resolve_path(v3, workdir) for k, v3 in dict(v2).items()
            }
    return reads


def preprocess_config(config, workdir):
    if "_REF" not in config:
        config["_REF"] = _preprocess_reference(config, workdir)
    if "_READS" not in config:
        config["_READS"] = _preprocess_reads(config, workdir)
