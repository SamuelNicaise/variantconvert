# -*- coding: utf-8 -*-
"""
@Author: Samuel Nicaise
@Version: v1.2.2
"""
import appdirs
import os

__version__ = "1.2.2"
__user_config__ = appdirs.AppDirs("variantconvert", version=__version__).user_config_dir
__default_config__ = os.path.join(os.path.dirname(__file__), "configs")
