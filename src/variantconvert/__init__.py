# -*- coding: utf-8 -*-
"""
@Author: Samuel Nicaise
@Version: v1.2.2
"""
import appdirs
import os

__version__ = "1.2.2"

__default_config__ = os.path.join(os.path.dirname(__file__), "configs")

__default_local_config__ = appdirs.AppDirs("variantconvert", version=__version__).user_config_dir

if os.path.isdir(os.environ.get("VARIANTCONVERT_CONFIG", "")):
    __local_config__ = os.environ["VARIANTCONVERT_CONFIG"]
else:
    __local_config__ = __default_local_config__
