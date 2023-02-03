# -*- coding: utf-8 -*-
"""
"""

import logging as log
import os
import pytest
import sys

sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
from variantconvert.config import get_nested_dic


def test_get_nested_dic():
    dic = {}
    target_key = "COLUMNS.INFO.ID.Test"
    expected = {"COLUMNS": {"INFO": {"ID": {"Test": ""}}}}
    res1 = get_nested_dic(dic, target_key)
    assert res1 == expected

    target_key = "COLUMNS.FILTER"
    default_val = "Pass"
    # fmt: off
    expected = {'COLUMNS': {
                        'INFO': {
                            'ID': {
                                'Test': ''
                            }
                        },
                        'FILTER': 'Pass'
                    } 
                }
    # fmt: on
    res2 = get_nested_dic(res1, target_key, default_value=default_val)
    assert res2 == expected


def test_change_config(tmp_path):
    pass
