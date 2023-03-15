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

    def manage_alt(self, row,ref):
        alt=""
        self.del_variant = False

        for i in self.sample_alt:
            values_sample = self.sample_alt.get(i)
            all_1 = values_sample.get(row)[0]
            all_2 = values_sample.get(row)[1]

            if all_1 == "-" or all_2 == "-" :
                if self.del_variant == False :
                    #je comprend pas pk ça rentre dedans
                    self.del_variant = True
                    ref=self.nuc_del+ref
                    #Changer le alt qui existe
                    if alt != "" :
                        alt_nw = ""
                        if alt == ".":
                            alt = ""
                        for str_alt in alt:

                            if str_alt != ',':
                                alt_nw = alt_nw + self.nuc_del + str_alt
                            else:
                                alt_nw = alt_nw + str_alt

                        alt = alt_nw

            if self.del_variant == True :
                if all_1 == "-" :
                    all_1=""
                if all_2 == "-":
                    all_2=""
                all_1=self.nuc_del+all_1
                all_2=self.nuc_del+all_2

            alt=self.generate_alt(alt, all_1, all_2, ref)
        return alt

    def generate_alt(self, alt, all_1, all_2, ref):
        # print(self.nuc_del," ", all_1," ", all_2," ", ref, " ",self.del_variant, alt)

        if alt == "":
            # Loresque le ALT est vide
            if all_1 == ref and all_2 == ref and alt.find(all_1) == -1 and alt.find(all_2) == -1:
                alt = "."

            if all_1 != ref and alt.find(all_1) == -1:
                # They are no one match whith alt and REF
                alt = all_1

            if all_2 != ref and alt.find(all_2) == -1:
                # They are no one match whith alt and REF
                alt = all_2

        elif alt != "":
            # Lorsque le ALT n'est pas vide
            # print("allele 1 ", all_1, " l'alleles 2 ", all_2, "la reff ", ref, "# le alt ", alt,
            #       " dans la boucle!!!!!!!!!!!!!!!!")

            if all_1 != ref and alt.find(all_1) == -1:
                # They are no one match whith alt and REF
                if alt.find(".") == 0:
                    alt = all_1
                else:
                    alt = alt + "," + all_1

            if all_2 != ref and alt.find(all_2) == -1:
                # They are no one match whith alt and REF
                if alt.find(".") == 0:
                    alt = all_2
                else :
                    alt = alt + "," + all_2

        if self.del_variant == True and len(all_1) == 1 and alt != "":
            if not self.search_del(alt) :
                alt = alt+","+all_1

        elif self.del_variant == True and len(all_2) == 1 and alt != "":
            if not self.search_del(alt) :
                alt = alt+","+all_2

        #gestion des deletions pour les trouvers

        return alt

    def search_del(self, alt):
        split_alt=alt.split(",")
        for i in split_alt :
            if len(i) == 1:
                return True
        #on split en plusieur le alt et juste regarder la longueur de chaque element pour voir si la délétion a était ajouté
        return False


    def _define_gt(self,row,ref):
        """If the nucleotide from the haplotype are find in REF we have a 0 and if it's find in alt it's 1. "/" correspond unphased genotype"""
        #Attention pour créer les reff on ne se base pas sur la méthode _get_alt_check_reff
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
            
            data = self.snp_data.astype(str).to_dict()

            error_reff=0
            for i in range(self.snp_data.shape[0]):
                # Ignorer line et on fait une fonction va afficher à la fin le variant voir pour faire un to string qui affiche toute les infos du variant à la fin
                var = Variant()
                self.nuc_del=""


                for vcf_col in ["#CHROM", "POS", "ID", "REF", "QUAL"]:
                    #Ici on va recuperer les infos de chaque colone 
                    col = self.config["VCF_COLUMNS"][vcf_col]

                    if is_helper_func(col):
                        #Si dans la config à coté du nom de colone on a mis en clé helper fonction ça va alors embrigader tout le reste. Ici ce sera pour notre ref
                        # col[1] is a function name, col[2] its list of args
                        # the function named in col[1] has to be callable from this module
                        func = helper.get(col[1]) #transformation du string en fonction python
                        args = [data[c][i] for c in col[2:]]
                        result = func(*args)
                        self.nuc_del = helper.nuc_del.upper()
                         #* permet de deballer tous les arguments qui sont enregistrer dans cette variable 
                        if result == "empty" :
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
                var.alt = self.manage_alt(i, var.ref)
                var.samples_gt = self._define_gt(i,var.ref)
                
                if var.ref == "":
                    continue

                if self.del_variant == True :
                    var.ref=self.nuc_del+var.ref

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

                if var.ref == "":
                    raise ValueError("They are no reff for variant", var.chrom, " ", var.pos)

                vcf.write("\t".join(line) + "\n")

        print(helper.error_value)
        #C'est pour connaitre le nombre de reff qui ne sont pas trouver dans le VCF

      

if __name__ == '__main__':
    convert = VcfFromSnp("/home1/BAS/hameaue/variant_convert/variantconvert/configs/config_snp.json")
    convert.convert("/home1/BAS/hameaue/TEST02_habibd/file_BBS5/Full_data_clean.csv", "/home1/BAS/hameaue/TEST02_habibd/file_BBS5/vcf_unphased.vcf")


