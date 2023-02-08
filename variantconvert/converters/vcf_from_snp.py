import pandas as pd
import os
import sys

from abstract_converter import AbstractConverter
sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
from helper_functions import HelperFunctions
from commons import create_vcf_header

class VcfFromSnp(AbstractConverter) :

    def data(self):
        ref_data=pd.read_csv('dbsnp.hg19.vcf',sep= '\t', index_col=0)

    def convert(self, input_path, output_path):
        print("hey")
        

if __name__ == '__main__':
    convert = VcfFromSnp("/home1/BAS/hameaue/variant_convert/variantconvert/configs/config_snp.json")
    convert.convert("/home1/BAS/hameaue/TEST02_habibd/cut_data.csv", "/home1/BAS/hameaue/TEST02_habibd/cut_data_unphased.vcf")


