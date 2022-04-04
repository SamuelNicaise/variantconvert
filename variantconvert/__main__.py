# -*- coding: utf-8 -*-
"""
@Goal: Expand Celine Besnard's script with infinite conversion abilities between vcf and various other formats
@Author: Samuel Nicaise
@Date: 23/11/2021

Prerequisites: pandas, pyfasta (https://github.com/brentp/pyfasta)

Usage examples:
#TSV (Decon) to VCF
python /home1/L/NGS/BIO_INFO/BIO_INFO_Sam/scripts/tsvConversion/file_converter.py -i /home1/BAS/nicaises/Tests/deconconverter/200514_NB551027_0724_AHTWHHAFXY.DECON_results_all.txt -o vcf_from_decon.vcf -fi tsv -fo vcf -c /home1/L/NGS/BIO_INFO/BIO_INFO_Sam/scripts/tsvConversion/fileconversion/config_decon.json

#AnnotSV3 to VCF
python /home1/L/NGS/BIO_INFO/BIO_INFO_Sam/scripts/tsvConversion/file_converter.py -i /home1/BAS/nicaises/Tests/deconconverter/DECoN.20211207-183955_results_both.tsv -o /home1/BAS/nicaises/Tests/deconconverter/decon__annotsv3.vcf -fi annotsv -fo vcf -c /home1/L/NGS/BIO_INFO/BIO_INFO_Sam/scripts/tsvConversion/fileconversion/config_annotsv3.json

#Canoes BED to VCF
python /home1/L/NGS/BIO_INFO/BIO_INFO_Sam/scripts/tsvConversion/file_converter.py -i canoes_calling.bed -o vcf_from_canoes_bed.vcf -fi tsv -fo vcf -c /home1/L/NGS/BIO_INFO/BIO_INFO_Sam/scripts/tsvConversion/fileconversion/config_canoes_bed.json

----------------------------
Configfile guidelines (JSON)
----------------------------
1)[GENERAL] has 3 important fields
	#source format name: will show up in VCF meta fields
	#skip_rows: how many rows to skip before we reach indexes.
	This script cannot handle a tsv with unnamed columns (beds are fine)
	#unique_variant_id: useful in multisample files. List the
	columns that are needed to uniquely identify a variant.
2) [VCF_COLUMNS] describe the columns that will go in your VCF
	key: column name in VCF ; value: column name in source format
3) [COLUMNS_DESCRIPTION] describe the tsv columns
	Type and Description fields will be used in the VCF header
4) read HelperFunctions docstring


#TODO: add argument mode to change config files (particularly genome)
#TODO: add argument mode to deal with an entire folder of varank files (or varank files in general)
#TODO: add a Variant class and a VCF class instead of writing files line by line
"""
from __future__ import division
from __future__ import print_function

import argparse
import os
import sys

from os.path import join as osj

from commons import set_log_level
from converter_factory import ConverterFactory


if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="python fileconverter.py")
    parser.add_argument("-i", "--inputFile", type=str, required=True, help="Input file")
    parser.add_argument(
        "-o", "--outputFile", type=str, required=True, help="Output file"
    )
    parser.add_argument(
        "-fi", "--inputFormat", type=str, required=True, help="Input file format"
    )
    parser.add_argument(
        "-fo", "--outputFormat", type=str, required=True, help="Output file format"
    )
    parser.add_argument(
        "-c",
        "--configFile",
        type=str,
        required=True,
        help="JSON config file describing columns. See script's docstring.",
    )
    parser.add_argument(
        "-v", "--verbosity", type=str, default="info", help="Verbosity level"
    )
    args = parser.parse_args()

    # varank_tester = type(
    #     "obj",
    #     (object,),
    #     {
    #         "inputFile": osj(
    #             os.path.dirname(__file__),
    #             "..",
    #             "..",
    #             "examples",
    #             "fam21_ACE2105946_BBS_allVariants.rankingByGene.tsv",
    #         ),
    #         "outputFile": osj(
    #             os.path.dirname(__file__), "..", "..", "examples", "varank_test.vcf"
    #         ),
    #         "inputFormat": "varank",
    #         "outputFormat": "vcf",
    #         "configFile": osj(
    #             os.path.dirname(__file__), "..", "configs", "config_varank.json"
    #         ),
    #         "verbosity": "debug",
    #     },
    # )
    # decon_tester = type(
    #     "obj",
    #     (object,),
    #     {
    #         "inputFile": osj(
    #             os.path.dirname(__file__),
    #             "..",
    #             "..",
    #             "examples",
    #             "DECON.20220329-083144.Design_results_all.txt",
    #         ),
    #         "outputFile": osj(
    #             os.path.dirname(__file__), "..", "..", "examples", "decon_test.vcf"
    #         ),
    #         "inputFormat": "tsv",
    #         "outputFormat": "vcf",
    #         "configFile": osj(
    #             os.path.dirname(__file__), "..", "configs", "config_decon.json"
    #         ),
    #         "verbosity": "debug",
    #     },
    # )
    # annotsv_tester = type(
    #     "obj",
    #     (object,),
    #     {
    #         "inputFile": osj(
    #             os.path.dirname(__file__),
    #             "..",
    #             "..",
    #             "examples",
    #             "DECoN.20211222-135425_results_all.AnnotSV.tsv",
    #         ),
    #         "outputFile": osj(
    #             os.path.dirname(__file__),
    #             "..",
    #             "..",
    #             "examples",
    #             "decon_annotsv_test.vcf",
    #         ),
    #         "inputFormat": "annotsv",
    #         "outputFormat": "vcf",
    #         "configFile": osj(
    #             os.path.dirname(__file__), "..", "configs", "config_annotsv3.json"
    #         ),
    #         "verbosity": "debug",
    #     },
    # )
    # bed_tester = type(
    #     "obj",
    #     (object,),
    #     {
    #         "inputFile": osj(
    #             os.path.dirname(__file__),
    #             "..",
    #             "..",
    #             "examples",
    #             "canoes.bed",
    #         ),
    #         "outputFile": osj(
    #             os.path.dirname(__file__), "..", "examples", "canoes_bed_test.vcf"
    #         ),
    #         "inputFormat": "tsv",
    #         "outputFormat": "vcf",
    #         "configFile": osj(
    #             os.path.dirname(__file__), "..", "configs", "config_canoes_bed.json"
    #         ),
    #         "verbosity": "debug",
    #     },
    # )
    # args = annotsv_tester

    set_log_level(args.verbosity)
    if args.inputFormat.lower() == "decon":
        print("[ERROR] DECON is handled as a TSV conversion. Use 'tsv' as input format")
        sys.exit()
    factory = ConverterFactory()
    converter = factory.get_converter(
        args.inputFormat.lower(), args.outputFormat.lower(), args.configFile
    )
    if args.inputFormat == "varank":
        converter.set_coord_conversion_file(
            osj(
                os.path.dirname(__file__),
                "..",
                "examples",
                "VCF_Coordinates_Conversion.tsv",
            )
        )
    converter.convert(args.inputFile, args.outputFile)
