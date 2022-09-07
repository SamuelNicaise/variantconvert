# -*- coding: utf-8 -*-

from __future__ import division
from __future__ import print_function

import logging as log
import os
import pandas as pd
import sys

from converters.abstract_converter import AbstractConverter

sys.path.append("..")
from commons import create_vcf_header, is_helper_func, clean_string
from helper_functions import HelperFunctions


class VcfFromBreakpoints(AbstractConverter):
    """Made for file formats such as the TSV outputs of STAR-Fusion and ARRIBA.
    Other converters (for now) are not able to generate a VCF containing breakpoints.
    
    Each input line will result in two VCF lines, one for each side of the breakpoint.
    """
    def _init_dataframe(self):
        self.df = pd.read_csv(
            self.filepath,
            skiprows=self.config["GENERAL"]["skip_rows"],
            sep="\t",
            low_memory=False,
        )
        self.df.reset_index(drop=True, inplace=True)
        self.df.fillna(".", inplace=True)
        log.debug(self.df)
        self.df["__!UNIQUE_VARIANT_ID!__"] = self.df.apply(
            lambda row: self._get_unique_variant_id(row), axis=1
        )
        log.debug(self.df)

    def _get_sample_list(self):
        # is the file multisample?
        if self.config["VCF_COLUMNS"]["SAMPLE"] != "":
            sample_list = []
            for sample in self.df[self.config["VCF_COLUMNS"]["SAMPLE"]].unique():
                sample_list.append(sample)
            return sample_list
        else:
            return [os.path.basename(self.output_path)]

    def _get_unique_variant_id(self, row):
        id = []
        for col in self.config["GENERAL"]["unique_variant_id"]:
            id.append(str(row[col]))
        return "_".join(id)

    def _get_unique_id_to_index_list(self, data):
        id_dic = {}
        for k, v in data["__!UNIQUE_VARIANT_ID!__"].items():
            if v not in id_dic:
                id_dic[v] = [k]
            else:
                id_dic[v].append(k)
        return id_dic

    def convert(self, tsv, output_path):
        log.info("Converting to vcf from annotSV using config: " + self.config_filepath)

        self.filepath = tsv
        self.output_path = output_path
        self._init_dataframe()
        sample_list = self._get_sample_list()
        helper = HelperFunctions(self.config)

        with open(output_path, "w") as vcf:
            vcf_header = create_vcf_header(tsv, self.config, sample_list, True)
            for l in vcf_header:
                vcf.write(l + "\n")

            data = self.df.astype(str).to_dict()
            # In some variant callers, output files contain a list of variant-sample associations
            # so the same variant can be on multiple lines
            # __!UNIQUE_VARIANT_ID!__ allows to identify such variants and only add them to the VCF once
            already_seen_variants = set()
            unique_id_to_index_list = self._get_unique_id_to_index_list(data)
            for i in range(len(data["__!UNIQUE_VARIANT_ID!__"])):
                if data["__!UNIQUE_VARIANT_ID!__"][i] in already_seen_variants:
                    continue

                lines = [[], []] # left side of the breakpoint, right side of the breakpoint

                for vcf_col in ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL"]:
                    col = self.config["VCF_COLUMNS"][vcf_col]

                    if vcf_col == "ID" and col == "":
                        #special override to name breakends
                        lines[0].append("bnd_" + str(i*2))
                        lines[1].append("bnd_" + str(i*2+1))
                        continue

                    if is_helper_func(col):
                        # col[1] is a function name, col[2] its list of args
                        # the function named in col[1] has to be callable from this module
                        func = helper.get(col[1])
                        args = [data[c][i] for c in col[2:]]
                        result = func(*args)
                        if len(result) != 2:
                            raise ValueError("HELPER_FUNCTIONS used with vcf_from_breakpoints.py are expected to return a tuple of len 2. Got instead:" + str(result))
                        lines[0].append(result[0])
                        lines[1].append(result[1])

                    elif col == "":
                        lines[0].append(".")
                        lines[1].append(".")
                    else:
                        lines[0].append(data[col][i])
                        lines[1].append(data[col][i])

                # Cutting-edge FILTER implementation
                lines[0].append("PASS")
                lines[1].append("PASS")

                left_info_field = ["SVTYPE=BND", "MATEID=" + "bnd_" + str(i*2+1)]
                right_info_field = ["SVTYPE=BND", "MATEID=" + "bnd_" + str(i*2)]
                for vcf_col, tsv_col in self.config["VCF_COLUMNS"]["INFO"].items():
                    if is_helper_func(tsv_col):
                        func = helper.get(tsv_col[1])
                        args = [data[c][i] for c in tsv_col[2:]]
                        s = vcf_col + "=" + func(*args)
                    else:
                        s = vcf_col + "=" + data[tsv_col][i]
                    left_info_field.append(clean_string(s))
                    right_info_field.append(clean_string(s))
                lines[0].append(";".join(left_info_field))
                lines[1].append(";".join(right_info_field))

                vcf_format_fields = []
                tsv_format_fields = []
                for vcf_col, tsv_col in self.config["VCF_COLUMNS"]["FORMAT"].items():
                    vcf_format_fields.append(vcf_col)
                    tsv_format_fields.append(tsv_col)
                lines[0].append(":".join(vcf_format_fields))
                lines[1].append(":".join(vcf_format_fields))

                # monosample input
                if len(sample_list) == 1:
                    sample_field = []
                    for key, val in self.config["VCF_COLUMNS"]["FORMAT"].items():
                        if key == "GT" and val == "":
                            sample_field.append("0/1")
                            continue
                        sample_field.append(data[val][index])
                    lines[0].append(":".join(sample_field))
                    lines[1].append(":".join(sample_field))
                # multisample input
                else:
                    sample_field_dic = {}
                    # If the variant exists in other lines in the source file, fetch their sample data now
                    for index in unique_id_to_index_list[
                        data["__!UNIQUE_VARIANT_ID!__"][i]
                    ]:
                        sample_field = []
                        for key, val in self.config["VCF_COLUMNS"]["FORMAT"].items():
                            if key == "GT" and val == "":
                                sample_field.append("0/1")
                                continue
                            sample_field.append(data[val][index])
                        sample_field_dic[
                            data[self.config["VCF_COLUMNS"]["SAMPLE"]][index]
                        ] = ":".join(sample_field)

                    for sample in sample_list:
                        if sample in sample_field_dic:
                            lines[0].append(sample_field_dic[sample])
                            lines[1].append(sample_field_dic[sample])
                        else:
                            if "GT" in self.config["VCF_COLUMNS"]["FORMAT"]:
                                if len(self.config["VCF_COLUMNS"]["FORMAT"]) == 1:
                                    # there's only GT. Avoid adding a trailing ":"
                                    empty = "./."
                                else:
                                    empty = "./.:" + ":".join(
                                        [
                                            "."
                                            for i in range(
                                                len(
                                                    self.config["VCF_COLUMNS"]["FORMAT"]
                                                )
                                                - 1
                                            )
                                        ]
                                    )
                            else:
                                empty = ":".join(
                                    [
                                        "."
                                        for i in range(
                                            len(self.config["VCF_COLUMNS"]["FORMAT"])
                                            - 1
                                        )
                                    ]
                                )
                            lines[0].append(empty)
                            lines[1].append(empty)
                    already_seen_variants.add(data["__!UNIQUE_VARIANT_ID!__"][i])
                
                #sort by chr/pos
                lines = sorted(lines, key=lambda x: (x[0], int(x[1])))
                for line in lines:
                    vcf.write("\t".join(line) + "\n")


if __name__ == "__main__":
    pass
