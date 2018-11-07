import numpy as np
import csv
import mysql.connector as connector
import os
import sys
import subprocess
from inparanoid_check import Inparanoid_check
from query import Query



# annotations currently in the local database:
# ['DB',
#  'SGDID',
#  'Symbol',
#  'Qualifier',
#  'GoId',
#  'GoReference',
#  'Evidence',
#  'WithOrFrom',
#  'Aspect',
#  'GeneProductName',
#  'Synonym',
#  'Type',
#  'taxon',
#  'Date',
#  'AnnotSource',
#  'PubMedID',
#  'Citation',
#  'GeneName',
#  'Feature',
#  'LiteratureTopic',
#  'SGDID',
#  'Ontology',
#  'AssociatedGenes',
#  'GOID',
#  'GoTerm',
#  'GoAspect',
#  'GoTermDef',
#  'ORF',
#  'SGDID',
#  'GENE_NAME',
# ]
# setting parameters to create the inparanoid_check object
my_local_ortho_database = "SGD1"
my_user = "root"
my_password = "CAmysql2280225964"
my_spec_name = "LP"
my_genome_file_name = "not_given"
my_gene_list = ["72701", "169754", "105970"]
my_annotation_list = ["GoId", "PubMedID"]
my_annot_dabase = "SGD1"

Inp_Check_obj = Inparanoid_check(my_spec_name, my_user, my_password, my_local_ortho_database, my_genome_file_name)

if not Inp_Check_obj.check_if_exists():
	Inp_Check_obj.run_inparanoid()
	parsed_inp_out = Inp_Check_obj.Inparanoid_Parser_()
	Inp_Check_obj.Build_Orthologs_database(parsed_inp_out)

my_query_obj = Query(my_spec_name, my_gene_list, my_annotation_list, my_annot_dabase, my_user, my_password )
my_query_obj.set_all_annot()
my_query_obj.set_all_genes()
my_query_obj.set_annot_dict()
my_query_obj.set_orthopairs()
my_query_obj.set_SGDID()
my_query_res=[]


for i in range(len(my_annotation_list)):
	my_query_res.append(my_query_obj.encode_to_string(my_query_obj.do_query()[i]))
	with open('results_'+ my_spec_name + '_' +str(i) + '.txt', 'w') as fp:
		fp.write('\n'.join('%s %s' % x for x in my_query_res[i]))




