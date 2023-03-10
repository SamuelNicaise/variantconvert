# variantconvert

The variantconvert is an extendable command-line tool for converting in VCF different file formats used to store genetic variant data. Currently, the following conversions are supported : 

- [AnnotSV](https://lbgi.fr/AnnotSV/) > VCF
- [STAR-Fusion](https://github.com/STAR-Fusion/STAR-Fusion) > VCF
- [Arriba](https://github.com/suhrig/arriba/) > VCF
- [DECoN](https://github.com/RahmanTeam/DECoN) > VCF
- BED (including [CANOES](https://github.com/bioinfo-chru-strasbourg/STARK-modules/tree/master/services/structuralvariation/canoes)) > VCF
- [VaRank](https://www.lbgi.fr/VaRank/) (single file or entire folder with varankBatch mode) > VCF

The project is still being developed and maintained. If you encounter any bug or wish to see new features added, please [open an issue](https://github.com/SamuelNicaise/variantconvert/issues).

# Installation

### Requirements

A minimal Python 3.8 installation is required, as well as the natsort, panda and pyfaidx Python modules

### Install
1) Setup an environment with Python >= 3.8. 
You can use the provided Conda [variantconvert.yml](https://github.com/SamuelNicaise/variantconvert/blob/master/variantconvert.yml) or the Python in the [Dockerfile](https://github.com/SamuelNicaise/variantconvert/blob/master/Dockerfile)
2) Do the following commands:
```
git clone https://github.com/SamuelNicaise/variantconvert.git
cd variantconvert
pip install -e .
```
3) Change the GENOME["path"] variable in configs/*.json to fit your local system. 
Indeed, some converters will require access to a valid reference genome in fasta format (used to fill in the REF column in cases where we only have the position without the reference base).
Alternatively, you can create your own config files in another folder.

# Usage
```
variantconvert --help 
```
Or if you did not use the `pip install` command above:
```
python variantconvert/__main__.py --help
```


____ 

# Documentation for AnnotSV user

<details> 
  <summary>Click to read documentation</summary>
  
### Creation of a VCF output file format with AnnotSV
To convert the output format from tsv to VCF, AnnotSV relies on the variantconvert tool. 

The variantconvert module distributed with AnnotSV can be used by setting the `-vcf` option to 1 in the AnnotSV command line.

### Requirements in the AnnotSV command line:
Different AnnotSV options are required to access to a VCF output:
-	From a "BED" or a "VCF" SV input file:
	- The user needs to define the `-SVinputInfo` option to 1 (to report in the tsv output file the 'ID', 'QUAL', 'FILTER'... fields).
-	From a "BED" SV input file:
	- The user needs to define the `-svtBEDcol` option (to report the SV type)
	- The `-samplesidBEDcol` option is highly recommended to use (else, the sample colum will be named "NA" (Non Attributed))  

### Method
Each SV from an AnnotSV tsv file is represented with 2 types of lines:
- An annotation on the “full” length of the SV. Every SV are reported, even those not covering a gene. 
- An annotation of the SV “split” by gene. This type of annotation gives an opportunity to focus on each gene overlapped by the SV. Thus, when a SV spans over several genes, the output will contain as many annotations lines as genes covered.

In the converted VCF, each SV is represented with only 1 line. All the annotations (full & split) are reported in the INFO field.
For one SV, all values from a same tsv output column are merged with a "|".

Example of a duplication overlapping 1 gene (1 full line + 1 split line in the tsv). The tsv output columns are represented in the INFO field in this way: 
```
AnnotSV_ID=21_35722427_35906593_DUP_1|21_35722427_35906593_DUP_1;SV_chrom=21|21;SV_start=35722427|35722427;SV_end=35906593|35906593;SV_lengt
h=184166|184166;SV_type=DUP|DUP;Annotation_mode=full&split;CytoBand=q22.12|q22.12;Gene_name=PPP1R2P2|PPP1R2P2;...
```
Warning: Space are replaced with an "_" in the output VCF

### GT warning:
It is to notice that if the GT is not given in input, the GT is set to “1/.” (using the variantconvert distributed by AnnotSV) or "0/1" (using the github variantconvert) for each SV in the VCF output file. Indeed, the considered SV has been called on at least one allele, but we don’t know the status of the second allele.
The user can configure this GT value directly in the variantconvert config files.

</detailS>


____

# Section for developers

<details> 
  <summary>Click to read documentation</summary>

## Adding new conversion formats

An intended goal of the project is to make it easy to add new formats to the conversion possibilities. 

Each conversion is described by a JSON config file with the following sections: 
- [GENERAL]
	- skip_rows: how many rows to skip before column indexes
	- unique_variant_id: A list of columns that are needed to uniquely identify a variant. Important for input files where a same variant can be on multiple lines. 

- [VCF_COLUMNS] maps input TSV columns to their corresponding VCF fields. 
  - Add or remove INFO fields at will to customize your output
  - When the equivalence is more complex than 1 input column = 1 VCF field ; you can create advanced HELPER_FUNCTION (explained below).
  
- [COLUMNS_DESCRIPTION]
  - Describes the input tsv columns to write the output VCF header. Column types can be inferred but it is usually safer to define them.

## HELPER_FUNCTION

They're defined in variantconvert/helperfunctions.py and called in your converter's config .json file. 

### To call a HELPER_FUNCTION

Use the following syntax in your .json: 
```
<vcf_field>: ["HELPER_FUNCTION", <function_name>, <tsv column 1>, <tsv column 2>...] # where tsv columns are the TSV fields sent as function input 
```

### To define a HELPER_FUNCTION

1. In HelperFunctions.__init__() , add *<function_name>* to the self.dispatcher dictionary
2. Add a new method in HelperFunctions class named as *<function_name>*, taking as parameters *<tsv column 1>, <tsv column 2>*... in the same order. Then you can use the full power of Python to do any data transformation you wish.

## If customizing a config file is not enough

variantconvert relies on Converter classes that are called by a ConverterFactory depending on the --inputFormat and --outputFormat parameters. 

You can create new Converter classes that will apply different transformations than the existing ones in variantconvert/converters/

They should inherit from the AbstractConverter class and be listed in the ConverterFactory class. That will make them automatically accessible from the command line. 
</details>
