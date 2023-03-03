# -*- coding: utf-8 -*-
"""
Use example for testing
variantconvert config --set GENOME.assembly=swag GENOME.path=/home1/DB/STARK/genomes/current/hg19.fa GENERAL.origin=myapp
"""

import itertools
import random
import logging as log
import os
import polars as pl
import secrets
import sys


from os.path import join as osj

sys.path.append(os.path.join(os.path.dirname(__file__), "."))
from commons import set_log_level


def load_vcf(vcf_path):
    VCF_HEADER = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
    header_metadata = []
    header_columns = ""
    content = []
    formats = []
    genotypes = []
    with open(vcf_path, "r") as f:
        for l in f:
            if l.startswith("##"):
                header_metadata.append(l)
            elif l.startswith("#"):
                header_columns = l
                if not header_columns.startswith(VCF_HEADER):
                    raise ValueError(
                        f"Expected VCF columns: {VCF_HEADER} followed by samples. Got instead: {header_columns}"
                    )
            elif len(l) > 0:
                l = l.rstrip().split("\t")
                content.append(l[:8])  # all columns until FORMAT (included)
                formats.append(l[8])  # FORMAT column
                genotypes.append(l[9:])  # all sample - variant associations
    return header_metadata, header_columns, content, formats, genotypes


def randomize_depths_in_sample_field(
    sample_field, format_field, random_min=10, random_max=300, seed=None
):
    """
    >>>randomize_depths_in_sample_field("0/1:30:20,10:yolo")
    "0/1:107:100,7:yolo"
    (using random)

    tl dr: randomize
    """
    if ":" in sample_field and ":" in format_field:
        sample_fields = sample_field.split(":")
        format_fields = format_field.split(":")
    elif ":" not in sample_field and ":" not in format_field:
        sample_fields = [sample_field]
        format_fields = [format_fields]
    else:
        raise ValueError(
            f"FORMAT and SAMPLE should have same number of fields. FORMAT: {format_field} ; SAMPLE: {sample_field}"
        )


def anonymize_depths(formats, genotypes, to_anonymize):
    if len(formats) != len(genotypes):
        raise ValueError(
            f"len(formats) = {len(formats)} ; len(genotypes) = {len(genotypes)} ; both should be equal"
        )

    # formats can be different between various lines: need to normalize columns to simplify computations
    genotype_list = []
    for l in range(len(genotypes)):
        for sample in genotypes[l]:
            genotype = {}
            formats_list = formats[l].split(":")
            for j in range(len(formats_list)):
                fields_list = sample.split(":")
                genotype[formats_list[j]] = fields_list[j]
            genotype_list.append(genotype)

    pl.Config.set_tbl_rows(20)
    df = pl.from_dicts(genotype_list)

    # Replace "." with null
    df = df.with_columns(
        pl.when(pl.col(pl.Utf8) == ".")
        .then(None)
        .otherwise(pl.col(pl.Utf8))  # keep original value
        .keep_name()
    )

    # df.select(pl.col('RE')) = df["RE"].cast(pl.Float64)
    # print(df)
    print(df.select(pl.max("RE")))
    df = df.with_columns(pl.col("RE").cast(pl.Float64))
    print(df)
    print(df.select(pl.max("RE")))
    # Note: this lambda doesn't do anything to null data, as wanted
    mean = df.select(pl.mean("RE"))[0, 0]
    std = df.select(pl.std("RE"))[0, 0]
    print(mean)
    print(std)
    df = df.with_columns(pl.col("RE").apply(lambda x: (x - mean) / (std)))
    print(df)


def randomize_genotypes(matrix, seed="default"):
    # verify that matrix is rectangular
    sublist_size = len(matrix[0])
    if sublist_size < 2:
        log.error(
            "Only one sample in input: genotypes can't be shuffled between samples. Output won't be truely anonymous."
        )
    for sublist in matrix:
        if len(sublist) != sublist_size:
            raise ValueError("Lines don't all have the same number of columns")

    # shuffle
    flattened_vcf = list(itertools.chain.from_iterable(matrix))
    if seed == "default":
        random.seed(
            secrets.token_hex()
        )  # ensure the seed is as random as possible; see Python doc for more info
    else:
        random.seed(seed)
    random.shuffle(flattened_vcf)

    # recreate matrix
    new_matrix = [
        flattened_vcf[x : x + sublist_size] for x in range(0, len(flattened_vcf), sublist_size)
    ]
    return new_matrix


def anonymize_samples(header_columns):
    """
    Takes "#CHROM..." vcf header line as input
    Returns it with all samples renamed as Sample1, Sample2... etc
    """
    cols = header_columns.split("\t")
    # columns were already checked in load_vcf() ; indexes can be used directly
    cols = cols[0:9] + [f"Sample{i+1}" for i in range(len(cols[9:]))]
    return "\t".join(cols)


def main_anonymize(args):
    set_log_level(args.verbosity)
    header_metadata, header_columns, content, formats, genotypes = load_vcf(args.input)
    header_columns = anonymize_samples(header_columns)
    anonymize_depths(formats, genotypes, ["RO, RR, RE"])
    genotypes = randomize_genotypes(genotypes)

    if len(content) != len(genotypes):
        raise ValueError(
            f"len(content) = {len(content)} ; len(genotypes) = {len(genotypes)} ; both should be equal"
        )

    with open(args.output, "w") as f:
        for l in header_metadata:
            f.write(l)
        f.write(header_columns + "\n")
        for i in range(len(content)):
            f.write("\t".join(content[i]) + "\t")
            f.write(formats[i] + "\t")
            f.write("\t".join(genotypes[i]) + "\n")


if __name__ == "__main__":
    print("hello")
    tester = type(
        "obj",
        (object,),
        {
            "input": osj(
                os.path.dirname(__file__),
                "..",
                "tests",
                "controls",
                "decon_annotsv_test.vcf",
            ),
            "output": osj(os.path.dirname(__file__), "..", "yolonaise.vcf"),
            "verbosity": "debug",
        },
    )

    main_anonymize(tester)
