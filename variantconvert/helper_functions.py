# -*- coding: utf-8 -*-

from __future__ import division
from __future__ import print_function

from commons import get_genome


class HelperFunctions:
    """
    For when you can't just convert columns by changing column names
    Steps needed:
    - define the helper function
    - update self.dispatcher so this class can redirect to the proper function
    - tell the config file you want to use a HELPER_FUNCTION with the following pattern:
            [HELPER_FUNCTION, <name of your function>, <arg1>, <arg2>, ..., <argN>]

    Example: I need a LENGTH value in my destination format,
    but my source file only has START and END columns.
    You would need:
    # somewhere in the module
            def get_length(start, end):
                    return str(end - start)
    #in this class __init__():
            self.dispatcher["get_length_from_special_format"]: get_length
    # in the JSON configfile
            LENGTH: [HELPER_FUNCTION, "get_length_from_special_format", START, END]
    """

    def __init__(self, config):
        self.config = config
        self.dispatcher = {
            "get_ref_from_decon": self.get_ref_from_decon,
            "get_alt_from_decon": self.get_alt_from_decon,
            "get_svlen_from_decon": self.get_svlen_from_decon,
            "get_info_from_annotsv": self.get_info_from_annotsv,
            "get_ref_from_canoes_bed": self.get_ref_from_canoes_bed,
            "get_alt_from_canoes_bed": self.get_alt_from_canoes_bed,
            "get_chr_from_breakpoint": self.get_chr_from_breakpoint,
            "get_pos_from_breakpoint": self.get_pos_from_breakpoint,
            "get_ref_from_breakpoint": self.get_ref_from_breakpoint,
            "get_alt_from_breakpoint": self.get_alt_from_breakpoint,
            "readable_starfusion_annots" : self.readable_starfusion_annots
        }

    def get(self, func_name):
        return self.dispatcher[func_name]

    def get_ref_from_decon(self, chr, start):
        f = get_genome(self.config["GENOME"]["path"])
        return f[chr][int(start) - 1].seq

    def get_ref_from_canoes_bed(self, chr, start):
        f = get_genome(self.config["GENOME"]["path"])
        print(self.config["GENOME"]["path"])
        print(f.keys())
        print("chr:", chr)
        return f["chr" + str(chr)][int(start) - 1].seq

    def get_ref_from_breakpoint(self, left_breakpoint, right_breakpoint):
        f = get_genome(self.config["GENOME"]["path"])

        left_chr = left_breakpoint.split(":")[0]
        if not left_chr.startswith("chr"):
            left_chr = "chr" + chr
        left_start = left_breakpoint.split(":")[1]

        right_chr = right_breakpoint.split(":")[0]
        if not right_chr.startswith("chr"):
            right_chr = "chr" + chr
        right_start = right_breakpoint.split(":")[1]

        return (f[left_chr][int(left_start) - 1].seq, f[right_chr][int(right_start) - 1].seq)

    def get_alt_from_breakpoint(self, left_breakpoint, right_breakpoint):
        left_ref, right_ref = self.get_ref_from_breakpoint(left_breakpoint, right_breakpoint)
        left_chr, left_pos, left_orientation = left_breakpoint.split(":")
        right_chr, right_pos, right_orientation = right_breakpoint.split(":")

        if left_orientation == "+":
            left_alt = f"{left_ref}[{right_chr}:{right_pos}["
        elif left_orientation == "-":
            left_alt = f"{left_ref}]{right_chr}:{right_pos}]"
        else:
            raise ValueError("Unexpected left_orientation:" + str(left_orientation))

        if right_orientation == "+":
            right_alt = f"]{left_chr}:{left_pos}]{right_ref}"
        elif right_orientation == "-":
            right_alt = f"[{left_chr}:{left_pos}[{right_ref}"
        else:
            raise ValueError("Unexpected right_orientation:" + str(right_orientation))

        return left_alt, right_alt

    @staticmethod
    def get_alt_from_decon(cnv_type_field):
        if cnv_type_field == "deletion":
            return "<DEL>"
        if cnv_type_field == "duplication":
            return "<DUP>"
        raise ValueError("Unexpected CNV.type value:" + str(cnv_type_field))

    @staticmethod
    def get_alt_from_canoes_bed(cnv_type_field):
        if cnv_type_field == "DEL":
            return "<DEL>"
        if cnv_type_field == "DUP":
            return "<DUP>"
        raise ValueError("Unexpected CNV.type value:" + str(cnv_type_field))

    @staticmethod
    def get_svlen_from_decon(start, end):
        return str(int(end) - int(start))

    @staticmethod
    def get_info_from_annotsv(info):
        """
        only used in attempts to convert annotsv files
        as if they were generic TSV (not recommended)
        """
        return "."

    @staticmethod
    def get_chr_from_breakpoint(left_breakpoint):
        return (left_breakpoint.split(":")[0], left_breakpoint.split(":")[0])

    @staticmethod
    def get_pos_from_breakpoint(left_breakpoint):
        return (left_breakpoint.split(":")[1], left_breakpoint.split(":")[1])

    @staticmethod
    def readable_starfusion_annots(annots):
        """
        input: '["Mitelman","ChimerKB","GUO2018CR_TCGA","DEEPEST2019","HGNC_GENEFAM","Cosmic","ChimerSeq","INTERCHROMOSOMAL[chr11--chr10]"]'
        output: 'Mitelman,ChimerKB,GUO2018CR_TCGA,DEEPEST2019,HGNC_GENEFAM,Cosmic,ChimerSeq,INTERCHROMOSOMAL[chr11--chr10]'
        """
        return ",".join([v[1:-1] for v in annots[1:-1].split(",")])
