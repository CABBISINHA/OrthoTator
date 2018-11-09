import numpy as np
import pandas as pd
import csv
import mysql.connector as connector
import os
import sys
import subprocess
from inparanoid_check import Inparanoid_check
from query import Query


def encode_to_string(encode_tuple_list): 
    aa_char_list = []
    for x in encode_tuple_list:
        aa_char_list.append([i.encode() for i in x])
    return(aa_char_list)

# setting parameters to create the inparanoid_check object
my_local_ortho_database = 'Inparanoid'
my_user_name = "root"
my_password = "CAmysql2280225964"
my_spec_name = "YARLI"
my_genome_file_name = "not_given"
my_gene_list = ["69905", "69598", "3304"]
# ["jgi|Rhoto_IFO0880_4|14897|fgenesh1_kg.6_#_1402_#_TRINITY_DN775_c0_g4_i1_mRNA", "jgi|Rhoto_IFO0880_4|12926|fgenesh1_pm.3_#_89_mRNA", "1502"]
my_annotation_list = ["GO_Id", "Phenotype"]
my_annot_dabase = "SGD"

# some LP genes
# "2478", "1777", "1502"
 
# some RT genes # RT doesn't work for now: parser need to be be fixed

# jgi|Rhoto_IFO0880_4|14897|fgenesh1_kg.6_#_1402_#_TRINITY_DN775_c0_g4_i1_mRNA
# jgi|Rhoto_IFO0880_4|12926|fgenesh1_pm.3_#_89_mRNA
# jgi|Rhoto_IFO0880_4|9502|gm1.5635_g_mRNA


# some YARLI genes
# ["69905", "69598", "3304"]


# set of available annotations:
# set(['Genetic_position', 'Start_coordinate', 'reference', 'DB', 'Feature', 'GO_Id', 'Chemical', 'Source',
# 'Feature_Name', 'Phenotype', 'Observed_condition', 'Allele', 'Parent_feature_name', 'Gene_name', 'Mutant_Type',
# 'Qualifier', 'Reporter', 'Type', 'Symbol', 'enzyme', 'Feature_name_Hit', 'Annotation_Source', 'Std_Gene_Name_Hit',
# 'Experiment_Type', 'Feature_name', 'PubMed_ID', 'Date', '2nd_SGDID', 'Feature_Name_Bait', 'Gene_Name', 'Feature_qualifier',
# 'Genetic_Physical_Interaction', 'With_or_From', 'Feature_Type', 'Synonym', 'Feature_type', 'Std_Gene_Name_Bait', 'Notes',
# 'Strand', 'GeneProductName', 'Taxon', 'Alias', 'GOID', 'SGDID', 'Gene', 'GO_Aspect', 'Sequence_version', 
# 'Standard_gene_name', 'literature_topic', 'Strain_Background', 'Reference', 'Stop_coordinate', 'DB_Reference', 'Citation', 
# 'EC_number', 'Evidence', 'Manually_curated', 'gene_name', 'Coordinate_version', 'Aspect', 'GO_Slim_term', 'ORF',
# 'Description', 'biochem_pathway', 'Chromosome', 'Details'])

IC = Inparanoid_check(my_spec_name, my_user_name, my_password, my_local_ortho_database, my_spec_name)

if not IC.check_if_exists():
	#Inp_Check_obj.run_inparanoid()
	parsed = IC.Inparanoid_Parser_()
	IC.Build_Orthologs_database(parsed)

Q = Query(my_spec_name, my_gene_list, my_annotation_list, my_local_ortho_database, my_annot_dabase, my_user_name, my_password)
my_query_res, my_query_res_conc = Q.do_query()


for i in range(len(my_query_res)):
    string_res = encode_to_string(my_query_res[i])
    with open('results_'+ my_spec_name + '_' + str(i) + '.txt', 'w') as file:
        file.writelines('\t'.join(x) + '\n' for x in string_res)

string_res = []
for i in range(len(my_query_res_conc)):
    string_res[i] = encode_to_string(my_query_res_conc[i])
    with open('results_conc_'+ my_spec_name + '_' + str(i) + '.txt', 'w') as file:
        file.writelines('\t'.join(x) + '\n' for x in string_res[i])

string_res = [x for x in string_res if len(x) > 0]
aa = []
for i in range(len(string_res)):
    aa.append(pd.DataFrame(string_res[i]))
aaa = reduce(lambda x, y: pd.merge(x, y, on = [0, 1, 2, 3]), aa)

writer = ExcelWriter('output_full.xlsx')
aaa.to_excel(writer,'Sheet1',index=False)
writer.save()







