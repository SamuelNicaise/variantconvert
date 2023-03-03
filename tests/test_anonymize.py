# -*- coding: utf-8 -*-
"""
"""

import logging as log
import os
import pytest
import sys

sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
from variantconvert.anonymize import anonymize_samples, anonymize_depths, randomize_genotypes


def test_anonymize_samples():
    header = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tbob01\tJEAN\tJ@cquel1ne"
    expected = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1\tSample2\tSample3"
    assert anonymize_samples(header) == expected
