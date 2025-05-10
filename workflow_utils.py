#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2025 Ye Chang yech1990@gmail.com
# Distributed under terms of the GNU license.
#
# Created: 2025-05-10 03:18


import os
from collections import defaultdict


def preprocess_config(config, workdir):
    def resolve_path(path):
        if os.path.isabs(path):
            return path
        else:
            path = os.path.expanduser(path)
            if os.path.isabs(path):
                return path
            else:
                return os.path.relpath(path, workdir)

    ref = config["reference"]
    for k, v in ref.items():
        ref[k] = [resolve_path(v2) for v2 in v]
    if "genome" not in ref:
        ref["genome"] = []
    # if contamination in ref but is empty, remove it
    if "contamination" in ref and len(ref["contamination"]) == 0:
        del ref["contamination"]

    reads = defaultdict(lambda: defaultdict(dict))
    for s, v in config[f"samples"].items():
        for i, v2 in enumerate(v["data"], 1):
            reads[str(s)][f"run{i}"] = {
                k: resolve_path(v3) for k, v3 in dict(v2).items()
            }

    return ref, reads
