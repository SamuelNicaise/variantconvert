import os
from os import path as op
from glob import glob
from setuptools import setup

from src import variantconvert

rel_default = op.relpath(variantconvert.__default_config__, start=op.dirname(__file__))
rel_user = op.relpath(variantconvert.__user_config__, start=op.dirname(__file__))

if not op.isdir(rel_user):
    os.makedirs(rel_user)

setup(
    data_files=[
        (
            rel_user,
            [f for f in glob(rel_default + "/**/*.json", recursive=True)],
        )
    ]
)
