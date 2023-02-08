import logging as log
import pandas as pd
import os
import sys
import re

from abstract_converter import AbstractConverter
sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
from helper_functions import HelperFunctions
from commons import create_vcf_header, is_helper_func
from variant import Variant

class VcfFromSnp(AbstractConverter) :

    def _init_dataframe(self):
        self.snp_data=pd.read_csv(
            self.filepath,
            sep= '\t',
            index_col=0
            )
        self.snp_data.reset_index(drop=True, inplace=True)

        

    def _get_sample_id(self):
        sample_list=[]
        for i in range(self.snp_data.shape[1]):
            string = self.snp_data.columns[i]
         
            found_allele = re.search('Top Alleles', string)

            if (found_allele != None):
                #keep top allele from raw and split them 
                split_allele = string.split('.')
                sample_list.append(split_allele[0])

        return sample_list

    def _add_values(self):
        #idée ajouter avec un générateur les valeurs qui sont nécéssaire
        pass
        


    def convert(self, input_path, output_path):
        log.debug("Converting to vcf from SNP using config: " + self.config_filepath)

        self.filepath=input_path
        self.output_path = output_path
        self._init_dataframe()

        sample_list = self._get_sample_id()
        helper = HelperFunctions(self.config)

        # # Permet de générer les collones
        with open(output_path, "w") as vcf:
            vcf_header = create_vcf_header(input_path, self.config, sample_list, True)
            for l in vcf_header:
                vcf.write(l + "\n")
            
            data = self.snp_data.astype(str).to_dict()

            for i in range(self.snp_data.shape[0]):
                line = [] 
                # Ignorer line et on fait une fonction va afficher à la fin le variant voir pour faire un to strig qui affiche toute les infos du variant à la fin
                var = Variant()

                for vcf_col in ["#CHROM", "POS", "ID", "REF", "QUAL"]:
                    col = self.config["VCF_COLUMNS"][vcf_col]

                    if is_helper_func(col):

                        # col[1] is a function name, col[2] its list of args
                        # the function named in col[1] has to be callable from this module
                        func = helper.get(col[1]) #transformation du string en fonction python
                        args = [data[c][i] for c in col[2:]]
                        result = func(*args) #* permet de deballer tous les arguments qui sont enregistrer dans cette variable 
                        if len(result) != 1:
                            raise ValueError(
                                "HELPER_FUNCTIONS used with vcf_from_snp.py are expected to return a tuple of len 1. Got instead:"
                                + str(result)
                            )
                        var.set_column(vcf_col, result)

                    elif col == "":
                        var.set_column(vcf_col, ".")
                    else:
                        var.set_column(vcf_col, data[col][i])




                    line = [
                        var.chrom,
                        var.pos,
                        # var.get_hash(),
                        # voir pour integrer id
                        var.ref.upper(),
                        #var.alt à faire
                        var.qual,
                    ]

                    line.append("PASS")
                print(line)





        
        
            


        

if __name__ == '__main__':
    convert = VcfFromSnp("/home1/BAS/hameaue/variant_convert/variantconvert/configs/config_snp.json")
    convert.convert("/home1/BAS/hameaue/TEST02_habibd/cut_data.csv", "/home1/BAS/hameaue/TEST02_habibd/cut_data_unphased.vcf")


