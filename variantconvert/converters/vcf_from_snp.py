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
from natsort import index_natsorted


class VcfFromSnp(AbstractConverter) :
    """Attention aux ALT vide, ils ne sont pas gérer pour le moment. C'est surement du à un probleme dans le is helper """
    def _init_dataframe(self):
        self.snp_data=pd.read_csv(
            self.filepath,
            sep= '\t',
            index_col=0
            )
        self.snp_data.reset_index(drop=True, inplace=True)
        

    def _get_sample_id(self):

        self.sample_list=[]
        self.sample_alt={}
        for i in range(self.snp_data.shape[1]):
            string = self.snp_data.columns[i]
         
            found_allele = re.search('Top Alleles', string)

            if (found_allele != None):
                row_alt={}
                """keep top allele from raw and split them """
                split_allele = string.split('.')
                self.sample_list.append(split_allele[0])

                for a in range(self.snp_data.shape[0]):
                    """For each id find we keep all haplotype """
                    row_alt[a]=self.snp_data.iloc[a,i]

                self.sample_alt[split_allele[0]]=row_alt

        return self.sample_list


    def _get_alt(self, row,ref):
        """Each haplotype in dico should create alt column, each nucleotide diffrents from REF are add to this string """
        alt=""
        for i in self.sample_alt:
            values_sample=self.sample_alt.get(i)

            all_1=values_sample.get(row)[0]
            all_2=values_sample.get(row)[1]

            if alt != "" : 
                if all_1 != ref and alt.find(all_1) == -1 :
                    #They are no one match whith alt and REF
                    alt = alt + "," + all_1
                if all_2 != ref and alt.find(all_2) == -1 :
                    #They are no one match whith alt and REF
                    alt = alt + "," + all_2
                    
            else :
                if all_1 != ref and alt.find(all_1) == -1 :
                    #They are no one match whith alt and REF (alt string are not empty)
                    alt = all_1
                if all_2 != ref and alt.find(all_2) == -1 :
                    #They are no one match whith alt and REF (alt string are not empty)
                    alt = all_2

            #TODO : verifier que ça prend bien en compte si on a pas de valeur est que c'est "-"
        return alt

    def _define_gt(self,row,ref):
        """If the nucleotide from the haplotype are find in REF we have a 0 and if it's find in alt it's 1. "/" correspond unphased genotype"""
        gt_samples={}

        for i in self.sample_alt:
            values_sample=self.sample_alt.get(i)
            gt=""

            for a in range(2) :
                allele = str(values_sample.get(row)[a])

                if allele == ref and gt == "" :
                    gt="0"
                
                elif allele != ref and gt == "" :
                    gt = "1"
                
                elif allele == ref and gt != "" :
                    gt = gt + "/0"
                
                elif allele != ref and gt != "" :
                    gt = gt + "/1"

            gt_samples[i] = gt
            list_gt=list(gt_samples.values())

        return list_gt

    def convert(self, input_path, output_path):
        log.debug("Converting to vcf from SNP using config: " + self.config_filepath)

        self.filepath=input_path
        self.output_path = output_path
        self._init_dataframe()

        sample_list = self._get_sample_id()
        helper = HelperFunctions(self.config)

        with open(output_path, "w") as vcf:
            vcf_header = create_vcf_header(input_path, self.config, sample_list, True)

            for l in vcf_header:
                vcf.write(l + "\n")
            
            # self.snp_data = self.snp_data.iloc[
            #     index_natsorted(self.snp_data[self.config["VCF_COLUMNS"]["#CHROM"]])
            # ]

            data = self.snp_data.astype(str).to_dict()

            error_reff=0
            for i in range(self.snp_data.shape[0]):
                # Ignorer line et on fait une fonction va afficher à la fin le variant voir pour faire un to string qui affiche toute les infos du variant à la fin
                var = Variant()

                for vcf_col in ["#CHROM", "POS", "ID", "REF", "QUAL"]:
                    col = self.config["VCF_COLUMNS"][vcf_col]

                    if is_helper_func(col):
                        # col[1] is a function name, col[2] its list of args
                        # the function named in col[1] has to be callable from this module
                        func = helper.get(col[1]) #transformation du string en fonction python
                        args = [data[c][i] for c in col[2:]]
                        result = func(*args) #* permet de deballer tous les arguments qui sont enregistrer dans cette variable 
                        if result == "empty" :
                            # print("They are not match with reff")
                            break

                        if len(result) != 1:
                            # print("func:", func)
                            # print("args:", args)
                            error_reff=error_reff+1
                            raise ValueError(
                                "HELPER_FUNCTIONS used with vcf_from_snp.py are expected to return a tuple of len 1. Got instead:"
                                + str(result)
                            )
                        # print(result)
                        var.set_column(vcf_col, result)

                    elif col == "":
                        var.set_column(vcf_col, ".")

                    else:
                        var.set_column(vcf_col, data[col][i])

                    var.ref=var.ref.upper()
                    var.alt=self._get_alt(i, var.ref)
                    var.samples_gt = self._define_gt(i,var.ref)

                    if var.alt == "" : 
                        continue
                    
                    if var.ref == "":
                        continue

                    line = [
                        "chr"+var.chrom,
                        var.pos,
                        var.id,
                        var.ref,
                        var.alt,
                        var.qual,
                        ]

                    line.append("PASS")
                    line.append(".")
                    line.append("GT")
                    line.extend(var.samples_gt)


                    vcf.write("\t".join(line) + "\n")

        print(helper.error_value)

      

if __name__ == '__main__':
    convert = VcfFromSnp("/home1/BAS/hameaue/variant_convert/variantconvert/configs/config_snp.json")
    convert.convert("/home1/BAS/hameaue/TEST02_habibd/Full_data_clean.csv", "/home1/BAS/hameaue/TEST02_habibd/vcf_unphased.vcf")


