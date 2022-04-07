# -*- coding: utf-8 -*-
"""
Test files are stored outside of the git project for confidentiality.
File structure for tests to work properly:

project_folder
    examples
        file_for_tests1.tsv
        file_for_tests2.tsv
    variantconvert
        .git
        tests
        config
        variantconvert
            __main__.py


For testing purposes, GENOME["path"] was changed in configs. Normally it is:
"path": "/home1/data/STARK/databases/genomes/current/hg19.fa"
"""
from __future__ import division
from __future__ import print_function

import logging as log
import os
import sys

from os.path import join as osj

sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
# from converter_factory import ConverterFactory
# from commons import set_log_level
from variantconvert.__main__ import main


def remove_if_exists(filepath):
    if os.path.exists(filepath):
        os.remove(filepath)


def test_varank_to_vcf():
    varank_tester = type(
        "obj",
        (object,),
        {
            "inputFile": osj(
                os.path.dirname(__file__),
                "..",
                "..",
                "examples",
                "fam21_SAMPLE_BBS_allVariants.rankingByGene.tsv",
            ),
            "outputFile": osj(
                os.path.dirname(__file__), "..", "..", "examples", "varank_test.vcf"
            ),
            "inputFormat": "varank",
            "outputFormat": "vcf",
            "configFile": osj(
                os.path.dirname(__file__), "..", "configs", "config_varank.json"
            ),
            "verbosity": "debug",
            "coordConversionFile": osj(
                os.path.dirname(__file__),
                "..",
                "..",
                "examples",
                "VCF_Coordinates_Conversion.tsv",
            ),
        },
    )
    remove_if_exists(varank_tester.outputFile)
    main(varank_tester)
    assert os.path.exists(varank_tester.outputFile)


def test_decon_to_vcf():
    decon_tester = type(
        "obj",
        (object,),
        {
            "inputFile": osj(
                os.path.dirname(__file__),
                "..",
                "..",
                "examples",
                "DECON.20220329-083144.Design_results_all.txt",
            ),
            "outputFile": osj(
                os.path.dirname(__file__), "..", "..", "examples", "decon_test.vcf"
            ),
            "inputFormat": "tsv",
            "outputFormat": "vcf",
            "configFile": osj(
                os.path.dirname(__file__), "..", "configs", "config_decon.json"
            ),
            "verbosity": "debug",
        },
    )
    remove_if_exists(decon_tester.outputFile)
    main(decon_tester)
    assert os.path.exists(decon_tester.outputFile)


def test_annotsv_to_vcf():
    annotsv_tester = type(
        "obj",
        (object,),
        {
            "inputFile": osj(
                os.path.dirname(__file__),
                "..",
                "..",
                "examples",
                "DECoN.20211222-135425_results_all.AnnotSV.tsv",
            ),
            "outputFile": osj(
                os.path.dirname(__file__),
                "..",
                "..",
                "examples",
                "decon_annotsv_test.vcf",
            ),
            "inputFormat": "annotsv",
            "outputFormat": "vcf",
            "configFile": osj(
                os.path.dirname(__file__), "..", "configs", "config_annotsv3.json"
            ),
            "verbosity": "debug",
        },
    )
    remove_if_exists(annotsv_tester.outputFile)
    main(annotsv_tester)
    assert os.path.exists(annotsv_tester.outputFile)


def test_bed_to_vcf():
    bed_tester = type(
        "obj",
        (object,),
        {
            "inputFile": osj(
                os.path.dirname(__file__),
                "..",
                "..",
                "examples",
                "canoes.bed",
            ),
            "outputFile": osj(
                os.path.dirname(__file__), "..", "..", "examples", "canoes_bed_test.vcf"
            ),
            "inputFormat": "tsv",
            "outputFormat": "vcf",
            "configFile": osj(
                os.path.dirname(__file__), "..", "configs", "config_canoes_bed.json"
            ),
            "verbosity": "debug",
        },
    )
    remove_if_exists(bed_tester.outputFile)
    main(bed_tester)
    assert os.path.exists(bed_tester.outputFile)


def test_celine_splice_data_to_vcf():
    celine_data_tester = type(
        "obj",
        (object,),
        {
            "inputFile": osj(
                os.path.dirname(__file__),
                "..",
                "..",
                "examples",
                "VUS_test.tsv",
            ),
            "outputFile": osj(
                os.path.dirname(__file__),
                "..",
                "..",
                "examples",
                "celine_splice_data.vcf",
            ),
            "inputFormat": "tsv",
            "outputFormat": "vcf",
            "configFile": osj(
                os.path.dirname(__file__),
                "..",
                "configs",
                "config_celine_splice_data.json",
            ),
            "verbosity": "debug",
        },
    )
    remove_if_exists(celine_data_tester.outputFile)
    main(celine_data_tester)
    assert os.path.exists(celine_data_tester.outputFile)


if __name__ == "__main__":
    test_varank_to_vcf()
