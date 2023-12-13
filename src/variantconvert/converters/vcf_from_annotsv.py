# -*- coding: utf-8 -*-

import logging as log
import os
import pandas as pd
import re
import sys
import time
from natsort import index_natsorted

from converters.abstract_converter import AbstractConverter

sys.path.append("..")
from commons import (
    create_vcf_header,
    info_string_to_dict,
    is_helper_func,
    is_int,
    is_float,
    remove_decimal_or_strip,
)

from helpers.annotsv_helper import AnnotSvHelper


class VcfFromAnnotsv(AbstractConverter):
    """
    Specificities compared to TSV:
    - vcf-like INFO field  in addition to other annotation columns
    - vcf-like FORMAT and <sample> fields
    - full/split annotations. Each variant can have one "full" and
    zero to many "split" annotations which result in additional lines in the file
    Very hard to deal with this with generic code --> it gets its own converter

    Sept 2022 update: add support for AnnotSV files obtained from bed files
    Those do not have REF, ALT, FORMAT and <sample_name> columns
    """

    def _build_input_dataframe(self):
        df = pd.read_csv(
            self.filepath,
            skiprows=self.config["GENERAL"]["skip_rows"],
            sep="\t",
            low_memory=False,
        )
        df.sort_values(
            [self.config["VCF_COLUMNS"]["#CHROM"], self.config["VCF_COLUMNS"]["INFO"]["SV_start"]],
            inplace=True,
        )
        df.reset_index(drop=True, inplace=True)

        sample_col = self.config["VCF_COLUMNS"]["SAMPLE"]
        if isinstance(sample_col, str) and sample_col != "":
            # avoid replacing "NA" sample by a dot
            df.loc[:, df.columns != sample_col] = df.loc[:, df.columns != sample_col].fillna(".")
            df.loc[:, sample_col] = df.loc[:, sample_col].fillna("NA")
        else:
            df.fillna(".", inplace=True)  # default empty value in VCF

        df = df.astype(str)
        log.debug(df)
        return df

    def _get_sample_list(self):
        if self.config["VCF_COLUMNS"]["FORMAT"] == "FORMAT":
            return self._get_sample_list_with_vcf_input()
        else:
            return self._get_sample_list_with_bed_input()

    def _get_sample_list_with_bed_input(self):
        samples_col = self.input_df[self.config["VCF_COLUMNS"]["SAMPLE"]]
        sample_list = []
        for cell in samples_col:
            if "," in cell:
                for sample in cell.split(","):
                    sample_list.append(sample)
            else:
                sample_list.append(cell)
        sample_list = list(set(sample_list))
        sample_list.sort()  # ensures output is always the same, despite using a set() above
        return sample_list

    def _get_sample_list_with_vcf_input(self):
        sample_list = []
        try:
            format_index = self.input_df.columns.to_list().index("FORMAT")
        except ValueError:
            raise ValueError(
                "ERROR: No FORMAT column detected. Variantconvert assumes you used a bed file as input for AnnotSV in that case. \nYou should use the 'annotsv3_from_bed.json' config, not 'annotsv3_from_vcf.json'"
            )

        for i in range(format_index + 1, self.input_df.columns.size):
            potential_sample_col = self.input_df[self.input_df.columns[i]]
            first_valid_index = potential_sample_col.first_valid_index()
            if self.is_sample_column(
                self.input_df.iloc[first_valid_index, i],
                self.input_df.iloc[first_valid_index, format_index],
            ):
                sample_list.append(self.input_df.columns[i])
            else:
                break

        log.debug(f"sample_list:{sample_list}")
        if not set(sample_list).issubset(self.input_df.columns):
            raise ValueError(
                f"When using an AnnotSV file generated from a VCF, all samples in {self.config['VCF_COLUMNS']['SAMPLE']} column are expected to have their own column in the input AnnotSV file"
            )
        return sample_list

    @staticmethod
    def is_sample_column(value, format):
        """
        Unperfect identification of sample columns. Assumed to be good enough within the context of AnnotSV files. Made to be used on columns directly after the FORMAT columns, knowing after sample columns should come the "Annotation" column (=full/split designation).
        Sample columns are arbitrarily identified because the "Samples" column
        """
        count = format.count(":")
        if count == 0:
            if format == "GT":
                # (.([|/].)+)  --> match if value is for example 0/1 or 1/0 or 1/2 or 0|1 or 0/1/2 etc...
                res = re.match(r"(.([|/].)+)", value)
                if res:
                    return True
                return False
            else:
                raise RuntimeError(
                    "Unable to determine which columns are sample columns.\nPlease verify that your input files are correct, then contact the developers at https://github.com/SamuelNicaise/variantconvert/issues"
                )
        else:
            if value.count(":") == count:
                return True
            return False

    def _build_input_annot_df(self):
        """
        remove FORMAT, Samples_ID, and each <sample> column
        TODO: remove vcf base cols ; INFO field
        """
        columns_to_drop = [v for v in self.sample_list]
        columns_to_drop += [v for v in self.main_vcf_cols]
        columns_to_drop.append(self.config["VCF_COLUMNS"]["SAMPLE"])
        columns_to_drop.append(self.config["VCF_COLUMNS"]["FORMAT"])

        if self.config["GENERAL"].get("keep_info", False) in ("true", "True"):
            log.debug(
                f"config['GENERAL']['keep_info'] was set to string:'{self.config['GENERAL']['keep_info']}' instead of a boolean. Converting value to True"
            )
            self.config["GENERAL"]["keep_info"] = True
        if self.config["GENERAL"].get("keep_info", False) == True:
            self.original_info_col = self.input_df[self.config["VCF_COLUMNS"]["INFO"]["INFO"]]
        columns_to_drop.append(self.config["VCF_COLUMNS"]["INFO"]["INFO"])
        for col in columns_to_drop:
            try:
                df = self.input_df.drop([col], axis=1)
            except KeyError:
                log.debug(f"Failed to drop column: {col}")

        # if self.config["GENERAL"].get("keep_info", False) == True:
        #     cols_without_info = df.columns.difference(["INFO"])
        #     df[cols_without_info] = df[cols_without_info].replace(
        #         ";", ",", regex=True
        #     )  # any ';' in annots will ruin the vcf INFO field
        # else:
        df = df.replace(";", ",", regex=True)  # any ';' in annots will ruin the vcf INFO field

        # TODO: check if CHROM col is in compliance with config ref genome (chrX or X)
        # if self.config["GENOME"]["vcf_header"][0].startswith("##contig=<ID=chr"):
        #     if not chrom.startswith

        if (
            self.config["VCF_COLUMNS"]["INFO"]["SVTYPE"] == ""
            or self.config["VCF_COLUMNS"]["INFO"]["SVTYPE"] not in df.columns
        ):
            raise ValueError(
                "SV_type column is required to turn an AnnotSV file into a VCF. Check if SV_type col is set in config or missing in your file.\n"
                + "If you generated your AnnotSV file from a bed, AnnotSV option -svtBEDcol is required."
            )
        return df

    def _merge_full_and_split(self, df):
        """
        input: df of a single annotSV_ID ; containing only annotations (no sample/FORMAT data)
        it can contain full and/or split annotations

        returns a single line dataframe with all annotations merged properly
        """
        if self.config["GENERAL"]["mode"] not in ("full&split", "full"):
            raise ValueError(
                "Unexpected value in json config['GENERAL']['mode']: "
                "only 'full&split' and 'full' mode are implemented."
            )
        # Do not keep 'base vcf col' in info field
        df = df.loc[
            :,
            [
                cols
                for cols in df.columns
                if cols not in ["ID", "REF", "ALT", "QUAL", "FILTER"] + self.sample_list
            ],
        ]

        annots = {}
        dfs = {}
        for typemode, df_type in df.groupby(self.config["VCF_COLUMNS"]["INFO"]["Annotation_mode"]):
            if typemode not in ("full", "split"):
                raise ValueError("Annotation type is assumed to be only 'full' or 'split'")
            dfs[typemode] = df_type

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
            # remove float decimal full rows
            for ann in dfs["full"].columns:
                annots[ann] = remove_decimal_or_strip(dfs["full"].loc[df.index[0], ann])

        if self.config["GENERAL"]["mode"] == "full":
            return annots

        # deal with split
        if "split" not in dfs.keys():
            return annots
        # each info field split
        for ann in dfs["split"].columns:
            transform = []
            if ann == self.config["VCF_COLUMNS"]["INFO"]["Annotation_mode"]:
                annots[ann] = self.config["GENERAL"]["mode"]
                continue
            # if you only keep full annotations split is lost advitam eternam
            # list of all values in each columns
            for splitval in dfs["split"][ann].tolist():
                transform.append(remove_decimal_or_strip(splitval))
            # we don't report split infos only if there are ALL equal to full row or they are stack to dot
            # if all([nq == annots[ann] for nq in transform]) or all(
            #     eq == "." for eq in transform
            # ):  # or ann in except_full_list:

            #    continue
            # In case of full and n split are different we keep values from all (more than 2 differencies)
            # else:
            values = [annots[ann]]
            values.extend(transform)
            annots[ann] = "|".join(values)
            # TODO in case of pipe already present in annotations change separator, maybe '+'

        # remove empty annots
        annots = {k: v for k, v in annots.items()}
        return annots

    def _build_info_dic(self):
        """
        Output: dictionary with key: annotsv_ID ; value: a key-value dictionary of all annotations
        This will be used to write the INFO field
        """
        input_annot_df = self._build_input_annot_df()
        annots_dic = {}
        id_col = self.config["VCF_COLUMNS"]["INFO"]["AnnotSV_ID"]
        supplemental_info_fields = []

        for variant_id, df_variant in input_annot_df.groupby(id_col):
            merged_annots = self._merge_full_and_split(df_variant)
            annots_dic[variant_id] = merged_annots

            if self.config["GENERAL"].get("keep_info", False) == True:
                info_dict = info_string_to_dict(
                    self.original_info_col.iloc[df_variant.index].iloc[0]
                )
                log.debug(f"info_dict before adding INFO: {info_dict}")
                for annot in info_dict:
                    if annot not in annots_dic[variant_id]:
                        annots_dic[variant_id][annot] = info_dict[annot]
                        if annot not in supplemental_info_fields:
                            supplemental_info_fields.append(annot)

                log.debug(f"info_dict after: {info_dict}")

        if self.config["GENERAL"].get("keep_info", False) == True:
            self.supplemental_header = self._build_supplemental_header(
                annots_dic, supplemental_info_fields
            )
        else:
            self.supplemental_header = []

        return annots_dic

    def _build_supplemental_header(self, annots_dic, additional_info_fields):
        """
        Create header strings for each value in additional_info_fields
        Types are inferred from the values in annots_dic

        Args:
            annots_dic (dict): has the following structure:
                                {
                                    variant1: {key1: val1, key2: val2},
                                    variant2: {key1: val1, key2: val2}
                                }
            additional_info_fields (set[str]): new info fields who need to be added to header.
                                            This list will be automatically completed with columns that
                                            are in annots_dic but not defined in config["COLUMNS_DESCRIPTION]

        Returns:
            list[str]: header strings, one for each additional_info_fields.
            Doesn't include newlines, as expected by commons.create_vcf_header(args)
        """
        supplemental_header = []
        known_descriptions = set(self.config["COLUMNS_DESCRIPTION"]["INFO"].keys())
        missing_annots = []

        # this function was originally made to rebuild the lost header of VCF annotations in "INFO" column in VCF>AnnotSV>VCF conversions.
        # at this point, we check if all INFO fields were defined in config.
        # if not, add a default header for them too.
        for variant, dic in annots_dic.items():
            for annot in dic:
                if annot not in known_descriptions and annot not in missing_annots:
                    missing_annots.append(annot)

        missing_annots = additional_info_fields + [
            v for v in missing_annots if v not in additional_info_fields
        ]
        log.debug(f"known desc:{missing_annots}")

        for field in missing_annots:
            # infer type
            if all([is_int(v.get(field, None)) for v in annots_dic.values()]):
                field_type = "Integer"
            elif all([is_float(v.get(field, None)) for v in annots_dic.values()]):
                field_type = "Float"
            else:
                field_type = "String"

            # infer number
            for info in annots_dic.values():
                if field in info:
                    if info[field] == None:
                        number = 0
                        field_type = "Flag"
                    else:
                        number = info[field].count(",") + 1

            if field in additional_info_fields:
                description = (
                    "Imported from the INFO field of the original VCF before AnnotSV annotation"
                )
            else:
                description = "Imported from AnnotSV"

            supplemental_header.append(
                "##INFO=<ID="
                + field
                + ",Number="
                + str(number)
                + ",Type="
                + field_type
                + ',Description="'
                + description
                + '">'
            )
        return supplemental_header

    # TODO: merge this with the other create_vcf_header method if possible
    # Making this method static is an attempt at making it possible to kick it out of the class
    @staticmethod
    def _create_vcf_header(input_path, config, sample_list, input_df, info_keys):
        header = []
        header.append("##fileformat=VCFv4.2")
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
        """
        Some columns (cols_to_init_now) are not directly linked to a Varank TSV column.
        They can be filled later, usually with a HELPER FUNC, but still have to be initiated in the data frame.
        """
        main_cols = []
        cols_to_init_now = []

        for vcf_col in ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER"]:
            tsv_col = self.config["VCF_COLUMNS"][vcf_col]
            if not isinstance(tsv_col, str):
                cols_to_init_now.append(vcf_col)
            elif tsv_col == "":
                cols_to_init_now.append(vcf_col)
            else:
                vcf_col = tsv_col
            main_cols.append(vcf_col)

        # adding all missing columns at once to avoid a PerformanceWarning
        self.input_df = pd.concat([self.input_df, pd.DataFrame(columns=cols_to_init_now)])
        self.input_df[cols_to_init_now] = ["." for i in range(len(cols_to_init_now))]

        return main_cols

    def _get_info_dict(self, df_variant, helper):
        # END, SVLEN and SVTYPE are reserved INFO keywords for SVs per VCF 4.3 specification
        # They are the only INFO field that have to be renamed
        ignored_cols = (
            self.main_vcf_cols + ["INFO", "SV_end", "SV_length", "SV_type"] + self.sample_list
        )

        for config_key, config_val in self.config["VCF_COLUMNS"]["INFO"].items():
            if is_helper_func(config_val):
                func = helper.get(config_val[1])
                args = [df_variant.get(c, None) for c in config_val[2:]]
                result = func(*args)
                df_variant[config_key] = result
        # TODO: helper func for the 3 info fields and done

        info_dic = {}
        for k, v in df_variant.items():
            if k not in ignored_cols:
                # if v == ".":
                #     info_dic[k] = None
                # else:
                info_dic[k] = remove_decimal_or_strip(v)
        return info_dic

    def convert(self, tsv, output_path):
        """
        Creates and fill the output file.

        For each annotSV_ID ; fetch all related lines of annotations in key value dics.
        A function identifies and merges the annotations.
        then we build a dictionary with 1 dictionary per annotSV_ID containing all the annotations
        This is then used to make the header and fill the INFO field.
        """
        log.debug("Converting to vcf from tsv using config: " + self.config_filepath)

        self.filepath = tsv
        helper = AnnotSvHelper(self.config)

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
            vcf_header = create_vcf_header(
                tsv,
                self.config,
                self.sample_list,
                supplemental_header=self.supplemental_header,
                fileformat="VCFv4.2",  # because no spaces in INFO fields for IGV support
            )
            for l in vcf_header:
                vcf.write(l + "\n")

            id_col = self.config["VCF_COLUMNS"]["INFO"]["AnnotSV_ID"]
            self.input_df = self.input_df.iloc[
                index_natsorted(self.input_df[self.config["VCF_COLUMNS"]["#CHROM"]])
            ]

            for df_variant in self.input_df.itertuples():
                # print(getattr(df_variant, "AnnotSV_ID"))
                df_variant = df_variant._asdict()  # for ease of dev

                if (
                    self.config["GENERAL"]["mode"] == "full"
                    and df_variant[self.config["VCF_COLUMNS"]["INFO"]["Annotation_mode"]] != "full"
                ):
                    continue

                # fill columns that need a helper func
                for config_key, config_val in self.config["VCF_COLUMNS"].items():
                    if config_key == "INFO":
                        continue  # dealt with by _get_info_dict()
                    elif config_key == "FILTER" and config_val == "":
                        df_variant[config_key] = "PASS"
                    elif is_helper_func(config_val):
                        func = helper.get(config_val[1])
                        args = [df_variant.get(c, None) for c in config_val[2:]]
                        result = func(*args)
                        df_variant[config_key] = result

                main_cols = "\t".join([df_variant[c] for c in self.main_vcf_cols])
                vcf.write(main_cols + "\t")

                info_list = []
                for k, v in self._get_info_dict(df_variant, helper).items():
                    if v != None:
                        info_list.append(k + "=" + v)
                    else:
                        info_list.append(k)  # deal with INFO flags
                vcf.write(";".join(info_list) + "\t")

                if self.config["VCF_COLUMNS"]["FORMAT"] != "":
                    sample_cols = "\t".join(
                        [df_variant[self.config["VCF_COLUMNS"]["FORMAT"]]]
                        + [df_variant[c] for c in self.sample_list]
                    )
                else:
                    sample_cols = "GT"
                    samples_with_variant = df_variant[self.config["VCF_COLUMNS"]["SAMPLE"]].split(
                        ","
                    )
                    for sample in self.sample_list:
                        if sample in samples_with_variant:
                            sample_cols += f"\t{self.config['GENERAL']['default_present_genotype']}"
                        else:
                            sample_cols += f"\t{self.config['GENERAL']['default_absent_genotype']}"
                vcf.write(sample_cols)
                vcf.write("\n")


# TODO: refonte de la gestion des full/split
# Si mode "full&split", au lieu d'ajouter des annotations, ajouter un variant par ligne avec la région correspondante
# TODO: deal with the particular case of INS where 1 full = 1 split --> regroup all annotations instead?
