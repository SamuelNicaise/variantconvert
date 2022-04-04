# -*- coding: utf-8 -*-

from __future__ import division
from __future__ import print_function

import logging as log
import os
import pandas as pd
import sys
import time

from converters.abstract_converter import AbstractConverter

sys.path.append("..")
from helper_functions import HelperFunctions


class VcfFromAnnotsv(AbstractConverter):
    """
    Specificities compared to TSV:
    - vcf-like INFO field  in addition to other annotation columns
    - vcf-like FORMAT and <sample> fields
    - full/split annotations. Each variant can have one "full" and
    zero to many "split" annotations which result in additional lines in the file
    Very hard to deal with this with generic code --> it gets its own converter
    """

    def _build_input_dataframe(self):
        df = pd.read_csv(
            self.filepath,
            skiprows=self.config["GENERAL"]["skip_rows"],
            sep="\t",
            low_memory=False,
        )
        df.sort_values(
            [self.config["VCF_COLUMNS"]["#CHROM"], self.config["VCF_COLUMNS"]["POS"]],
            inplace=True,
        )
        df.reset_index(drop=True, inplace=True)
        df.fillna(".", inplace=True)
        df = df.astype(str)
        log.debug(df)
        return df

    def _get_sample_list(self):
        samples_col = self.input_df[self.config["VCF_COLUMNS"]["SAMPLE"]]
        sample_list = []
        for cell in samples_col:
            if "," in cell:
                for sample in cell.split(","):
                    sample_list.append(sample)
            else:
                sample_list.append(cell)
        # print(samples_col)
        sample_list = list(set(sample_list))
        # print("sample_list:", sample_list)
        if not set(sample_list).issubset(self.input_df.columns):
            raise ValueError(
                "All samples in '" + samples_col + "' column are expected to "
                "have their own column in the input AnnotSV file"
            )
        return sample_list

    def _build_input_annot_df(self):
        """
        remove FORMAT, Samples_ID, and each <sample> column
        TODO: remove vcf base cols ; INFO field
        """
        columns_to_drop = [v for v in self.sample_list]
        columns_to_drop += [v for v in self.main_vcf_cols]
        columns_to_drop.append(self.config["VCF_COLUMNS"]["SAMPLE"])
        columns_to_drop.append(self.config["VCF_COLUMNS"]["FORMAT"])
        columns_to_drop.append(self.config["VCF_COLUMNS"]["INFO"]["INFO"])
        df = self.input_df.drop(columns_to_drop, axis=1)
        df = df.replace(
            ";", ",", regex=True
        )  # any ';' in annots will ruin the vcf INFO field
        return df

    def _merge_full_and_split(self, df):
        """
        input: df of a single annotSV_ID ; containing only annotations (no sample/FORMAT data)
        it can contain full and/or split annotations

        returns a single line dataframe with all annotations merged properly
        """
        if self.config["GENERAL"]["mode"] != "full&split":
            raise ValueError(
                "Unexpected value in json config['GENERAL']['mode']: "
                "only 'full&split' mode is implemented yet."
            )

        annots = {}
        dfs = {}
        for type, df_type in df.groupby(
            self.config["VCF_COLUMNS"]["INFO"]["Annotation_mode"]
        ):
            if type not in ("full", "split"):
                raise ValueError(
                    "Annotation type is assumed to be only 'full' or 'split'"
                )
            dfs[type] = df_type

        # deal with full
        if "full" not in dfs.keys():
            # still need to init columns
            if "split" not in dfs.keys():
                log.warning(
                    "Input does not include AnnotSV's 'Annotation_mode' column. This is necessary to know how to deal with annotations. The INFO field will be empty."
                )
                return {}
            for ann in dfs["split"].columns:
                annots[ann] = "."
        else:
            if len(dfs["full"].index) > 1:
                raise ValueError(
                    "Each variant is assumed to only have one single line of 'full' annotation"
                )
            for ann in dfs["full"].columns:
                annots[ann] = dfs["full"].loc[df.index[0], ann]

        # deal with split
        if "split" not in dfs.keys():
            return annots
        for ann in dfs["split"].columns:
            if ann == self.config["VCF_COLUMNS"]["INFO"]["Annotation_mode"]:
                annots[ann] = self.config["GENERAL"]["mode"]
                continue
            if annots[ann] != ".":
                continue  #'full' annot is always prioritized
            annots[ann] = ",".join(dfs["split"][ann].tolist())

        # remove empty annots
        annots = {k: v for k, v in annots.items() if v != "."}
        return annots

    def _build_info_dic(self):
        """
        Output: dictionary with key: annotsv_ID ; value: a key-value dictionary of all annotations
        This will be used to write the INFO field
        """
        input_annot_df = self._build_input_annot_df()
        # print(input_annot_df)
        annots_dic = {}
        id_col = self.config["VCF_COLUMNS"]["INFO"]["AnnotSV_ID"]
        for variant_id, df_variant in input_annot_df.groupby(id_col):
            merged_annots = self._merge_full_and_split(df_variant)
            annots_dic[variant_id] = merged_annots
        return annots_dic

    # TODO: merge this with the other create_vcf_header method if possible
    # Making this method static is an attempt at making it possible to kick it out of the class
    @staticmethod
    def _create_vcf_header(input_path, config, sample_list, input_df, info_keys):
        header = []
        header.append("##fileformat=VCFv4.3")
        header.append("##fileDate=%s" % time.strftime("%d/%m/%Y"))
        header.append("##source=" + config["GENERAL"]["origin"])
        header.append("##InputFile=%s" % os.path.abspath(input_path))

        if config["VCF_COLUMNS"]["FILTER"] in input_df.columns:
            for filter in set(input_df[config["VCF_COLUMNS"]["FILTER"]].to_list()):
                header.append("##FILTER=<ID=" + str(filter) + ',Description=".">')
        else:
            header.append('##FILTER=<ID=PASS,Description="Passed filter">')

        # identify existing values in header_dic...
        header_dic = {}
        header_dic["REF"] = set(input_df[config["VCF_COLUMNS"]["REF"]].to_list())
        header_dic["ALT"] = set(input_df[config["VCF_COLUMNS"]["ALT"]].to_list())
        # ... then for each of them, check if a description was given in the config
        for section in header_dic:
            for key in header_dic[section]:
                if key in config["COLUMNS_DESCRIPTION"][section]:
                    info_config = config["COLUMNS_DESCRIPTION"][section][key]
                    header.append(
                        "##"
                        + section
                        + "=<ID="
                        + key
                        + ',Description="'
                        + info_config["Description"]
                        + '">'
                    )
                else:
                    header.append(
                        "##"
                        + section
                        + "=<ID="
                        + key
                        + ',Description="Imported from '
                        + config["GENERAL"]["origin"]
                        + '">'
                    )

        # same as before, but for header elements that also have a type...
        header_dic = {}
        header_dic["INFO"] = info_keys
        header_dic["FORMAT"] = set()
        for format_field in input_df[config["VCF_COLUMNS"]["FORMAT"]].to_list():
            for format in format_field.split(":"):
                header_dic["FORMAT"].add(format)
        # ... then for each of them, check if a description was given in the config
        for section in header_dic:
            for key in header_dic[section]:
                if key in config["COLUMNS_DESCRIPTION"][section]:
                    info_config = config["COLUMNS_DESCRIPTION"][section][key]
                    header.append(
                        "##"
                        + section
                        + "=<ID="
                        + key
                        + ",Number=.,Type="
                        + info_config["Type"]
                        + ',Description="'
                        + info_config["Description"]
                        + '">'
                    )
                else:
                    header.append(
                        "##"
                        + section
                        + "=<ID="
                        + key
                        + ',Number=.,Type=String,Description="Imported from '
                        + config["GENERAL"]["origin"]
                        + '">'
                    )

        header += config["GENOME"]["vcf_header"]
        header.append(
            "\t".join(
                [
                    "#CHROM",
                    "POS",
                    "ID",
                    "REF",
                    "ALT",
                    "QUAL",
                    "FILTER",
                    "INFO",
                    "FORMAT",
                ]
                + sample_list
            )
        )
        return header

    def _get_main_vcf_cols(self):
        cols = [self.config["VCF_COLUMNS"]["#CHROM"]]
        cols.append(self.config["VCF_COLUMNS"]["POS"])
        cols.append(self.config["VCF_COLUMNS"]["ID"])
        cols.append(self.config["VCF_COLUMNS"]["REF"])
        cols.append(self.config["VCF_COLUMNS"]["ALT"])
        cols.append(self.config["VCF_COLUMNS"]["QUAL"])
        cols.append(self.config["VCF_COLUMNS"]["FILTER"])
        return cols

    def convert(self, tsv, output_path):
        """
        Creates and fill the output file.

        For each annotSV_ID ; fetch all related lines of annotations in key value dics.
        A function identifies and merges the annotations.
        then we build a dictionary with 1 dictionary per annotSV_ID containing all the annotations
        This is then used to make the header and fill the INFO field.

        Note: the "INFO" field from annotSV is discarded for now,
        because it only contains Decon annotations and they're useless.
        TODO: make an option to keep the "INFO" field in the annotations dictionary
        """
        log.info("Converting to vcf from tsv using config: " + self.config_filepath)

        self.filepath = tsv
        helper = HelperFunctions(self.config)

        self.input_df = self._build_input_dataframe()
        self.sample_list = self._get_sample_list()
        self.main_vcf_cols = self._get_main_vcf_cols()

        info_dic = self._build_info_dic()
        info_keys = set()
        for id, dic in info_dic.items():
            for k in dic:
                info_keys.add(k)

        # create the vcf
        with open(output_path, "w") as vcf:
            vcf_header = self._create_vcf_header(
                tsv, self.config, self.sample_list, self.input_df, info_keys
            )
            for l in vcf_header:
                vcf.write(l + "\n")

            id_col = self.config["VCF_COLUMNS"]["INFO"]["AnnotSV_ID"]
            for variant_id, df_variant in self.input_df.groupby(id_col):
                main_cols = "\t".join(df_variant[self.main_vcf_cols].iloc[0].to_list())
                vcf.write(main_cols + "\t")
                vcf.write(
                    ";".join([k + "=" + v for k, v in info_dic[variant_id].items()])
                    + "\t"
                )
                sample_cols = "\t".join(
                    df_variant[
                        [self.config["VCF_COLUMNS"]["FORMAT"]] + self.sample_list
                    ]
                    .iloc[0]
                    .to_list()
                )
                vcf.write(sample_cols + "\t")
                vcf.write("\n")


if __name__ == "__main__":
    pass
