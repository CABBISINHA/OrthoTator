import numpy as np
import pandas as pd
import csv
import mysql.connector as connector
import os
import sys
import subprocess
from Inparanoid_check import Inparanoid_check
from query import Query
from absl import flags
from absl import logging
from absl import app


FLAGS = flags.FLAGS

flags.DEFINE_string('query', None, 'name of the query genes file, this is a tab delimited file containing all query genes')
flags.DEFINE_string('annot', None, 'name of the annot file, this is a tab delimited file containing all annotations')
flags.DEFINE_string('species', None, 'Species name')
flags.DEFINE_string('genome', None, 'name of the FASTA formatted genome file')
flags.DEFINE_string('ortho_db', None, 'name of the local orthology databse')
flags.DEFINE_string('sgd_db', None, 'name of the local SGD databse')
flags.DEFINE_string('user', None, 'username to connect to local mySQL')
flags.DEFINE_string('password', None, 'password to connect to local mySQL')

def main(argv):
    IC = Inparanoid_check(FLAGS.species, FLAGS.user, FLAGS.password, FLAGS.ortho_db, FLAGS.genome)
    if not IC.check_if_exists():
    	IC.run_inparanoid()
    	parsed = IC.Inparanoid_Parser_()
    	IC.Build_Orthologs_database(parsed)

    with open(FLAGS.query,'r') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            all_genes = row

    with open(FLAGS.annot,'r') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            all_annots = row


    Q = Query(FLAGS.species, all_genes, all_annots, FLAGS.ortho_db, FLAGS.sgd_db, FLAGS.user, FLAGS.password)
    my_query_res, my_query_res_conc = Q.do_query()

    for i in range(len(my_query_res)):
        string_res = encode_to_string(my_query_res[i])
        with open('results_'+ FLAGS.species + '_' + str(i) + '.txt', 'w') as file:
            file.writelines('\t'.join(x) + '\n' for x in string_res)

    string_res = []
    for i in range(len(my_query_res_conc)):
        string_res.append( encode_to_string(my_query_res_conc[i]))
        with open('results_conc_'+ FLAGS.species + '_' + str(i) + '.txt', 'w') as file:
            file.writelines('\t'.join(x) + '\n' for x in string_res[i])

    string_res = [x for x in string_res if len(x) > 0]
    aa = []
    for i in range(len(string_res)):
        aa.append(pd.DataFrame(string_res[i]))
    aaa = reduce(lambda x, y: pd.merge(x, y, on = [0, 1, 2, 3]), aa)

    writer = pd.ExcelWriter('output_full.xlsx')
    aaa.to_excel(writer,'Sheet1',index=False)
    writer.save()



def encode_to_string(encode_tuple_list):
    aa_char_list = []
    for x in encode_tuple_list:
        aa_char_list.append([i.encode() for i in x])
    return(aa_char_list)


if __name__ == "__main__":
    app.run(main)
