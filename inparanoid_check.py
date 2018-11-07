import numpy as np
import csv
import mysql.connector as connector
import os
import sys
import subprocess

class Inparanoid_check(object):
    #TODO: Check if the input query species name is already in the database

    def __init__(self,spec_name, database_user, database_password, local_Inparanoid_database_name, genome_file_name = "not_given"):
        #######
        # Inputs:
        #genome_file_name: name of the FASTA formatted file suitable for inparanoid input
        #
        self.spec_name_ = spec_name
        self.genome_file_name_ = genome_file_name
        self.database_user_ = database_user
        self.database_password_ = database_password
        self.local_inparanoid_ = local_Inparanoid_database_name


    def check_if_exists(self):
        database_connection = connector.connect(host="localhost",
                 user = self.database_user_,
                 passwd = self.database_password_)
        cursor = database_connection.cursor()
        cursor.execute('USE ' + self.local_inparanoid_)

        cursor.execute("SHOW TABLES LIKE " + "'" + self.spec_name_ + "_SC" + "'")
        result = cursor.fetchall()
        cursor.close()

        if result == []:
            return False
        else:
            return True

    ##########
    # The following methods may be used in main if "check_if_exists" returns false
    ##########

    def run_inparanoid(self):
        #os.system("sh utils/run_inparanoid.sh")
        #TODO, first finilize the inparanoid bash script then complete it
        print("start running InParanoid")
        subprocess.check_call("./run_inparanoid.sh %s" %(str(self.genome_file_name_)), shell=True) # if not working consider passing shell= true
        print("done running InParanoid")


    def Inparanoid_Parser_(self):
        # This functions returns a python list that can be used for "Build_Orthologs_database"
        geneSpec1 = []
        geneSpec2 = []
        scoreSpec1 = []
        scoreSpec2 = []
        totalScore = []
        clusterID = 1

        with open('sqltable.' + self.genome_file_name_ + '-SC','r') as f:
            reader = csv.reader(f, delimiter='\t')
            ret = []
            for row in reader:
                if(np.int(row[0]) != clusterID):
                    for g1, s1 in zip(geneSpec1, scoreSpec1):
                        for g2, s2 in zip(geneSpec2, scoreSpec2):
                            currRow = [g1,g2,s1,s2,totalScore[-1]]
                            ret.append(currRow)

                    geneSpec1 = []
                    geneSpec2 = []
                    scoreSpec1 = []
                    scoreSpec2 = []
                    totalScore = []

                    clusterID += 1

                if (row[2] == self.genome_file_name_):
                    geneSpec1.append(row[4])
                    scoreSpec1.append(np.float(row[3]))
                    totalScore.append(np.int(row[1]))
                elif (row[2] == 'SC'):
                    geneSpec2.append(row[4])
                    scoreSpec2.append(np.float(row[3]))
                    totalScore.append(np.int(row[1]))

        return ret

        #TODO, convert R code written by Shayan into a python method

    def Build_Orthologs_database(self, python_list_parsed):
        database_connection = connector.connect(host="localhost",
                 user = self.database_user_,
                 passwd = self.database_password_)
        cursor = database_connection.cursor()
        cursor.execute('USE ' + self.local_inparanoid_)
        cursor.execute('CREATE TABLE ' + self.spec_name_ + '_SC (' +  self.spec_name_ +  '_genes VARCHAR(200), SC_genes VARCHAR(200), Score_g1 DOUBLE(2,2), Score_g2 DOUBLE(2,2), Score_in INT(6))')


        for row in python_list_parsed:
            cursor.execute("INSERT INTO " + self.spec_name_ + '_SC (' + self.spec_name_ + "_genes, SC_genes, Score_g1, Score_g2, Score_in) VALUES(%s, %s, %s, %s, %s)", row)

        database_connection.commit()
        cursor.close()

#############################################################################
