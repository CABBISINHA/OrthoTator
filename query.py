import csv
import mysql.connector as connector
import os
import sys
import subprocess

class Query:
    def __init__(self, species, gene_list, annot_list, inparanoid_dbName, SGD_dbName, db_user_name, db_password):
        self.spc = species
        self.glist = gene_list
        self.annlist = annot_list
        self.spc_table = species + '_' + 'SC'
        self.inparanoid_db = inparanoid_dbName
        self.SGD_db = SGD_dbName
        self.user_name = db_user_name
        self.password = db_password

        self.inparanoid_DB_connection = connector.connect(host="localhost",
                                                user=self.user_name,
                                                passwd=self.password,
                                                database=self.inparanoid_db)

        self.SGD_DB_connection = connector.connect(host="localhost",
                                                user=self.user_name,
                                                passwd=self.password,
                                                database=self.SGD_db)

        if not self.mapping_table_exists():
            self.build_SGDID_map_db()
        print "done_building_SGDID_Map"
        self.set_all_genes()
        print "done_setting_all_genes"
        self.set_ortho_map()
        print "done_setting_ortho_map"
        self.all_annots = self.set_annot_dict() # This also sets annot_dict
        print "done_setting_all_annots"


    def gene_list_to_char(self, list_of_genes):
        my_char="("
        for gene in list_of_genes:
            my_char = my_char + "'" + gene + "'" + ", "
        my_char = my_char[0:-2]
        my_char = my_char + ")"
        return my_char

    def set_all_genes(self):
        # this function gets the list of genes corresponding to the species in the orhology table
        cursor = self.inparanoid_DB_connection.cursor()
        cursor.execute("SELECT " + self.spc + '_genes' + " FROM " + self.spc_table + ";")
        gname_encode = cursor.fetchall()
        cursor.close()

        gname_str = []
        for gtuple in gname_encode:
            gname_str.append(gtuple[0].encode())
        self.all_genes = set(gname_str)


    def check_gene_exist(self, gene_name):
        exists = set(gene_name).issubset(self.all_genes)
        return exists

    def set_ortho_map(self):
        ## consider adding a threshold as input to filter orthologs if needed
        """
        This function retrieves the orthologs of the gene_list in SC.
        It reads the orthologs from a table in the database. the table is named: 'desired_species_name'_SC
        and the column is named 'desired_species_name'.

        It forms a dictonary, keys are genes in gene_list (those that have orthologs), value is another
        dictopnary with two keys (for now) 'ORF' and 'SGDID', value for each is an array containing all
        orthologs.
        """
        cursor = self.inparanoid_DB_connection.cursor()

        ortho_query = " WHERE " + self.spc + '_genes' + " IN " +  self.gene_list_to_char(self.glist) +";"
        cursor.execute("SELECT " + self.spc+ '_genes' + ", SC_genes FROM " + self.spc_table + ortho_query)
        my_ortho_pairs = cursor.fetchall()
        cursor.close()

        if my_ortho_pairs == []:
            print "None of the provided genes exist in the orthology table."
            return

        orthoMap = {}
        for p in my_ortho_pairs:
            key = p[0].encode()

            if key in orthoMap.keys():
                orthoMap[key]['ORF'].append(p[1].encode())
                orthoMap[key]['SGDID'].append(self.convert_to_SGDID(p[1].encode()))
            else:
                orthoMap[key] = {}
                orthoMap[key]['ORF'] = []
                orthoMap[key]['SGDID'] = []
                orthoMap[key]['ORF'].append(p[1].encode())
                orthoMap[key]['SGDID'].append(self.convert_to_SGDID(p[1].encode()))

        self.ortho_map = orthoMap

    def mapping_table_exists(self):
        cursor = self.SGD_DB_connection.cursor()

        cursor.execute("SHOW TABLES LIKE 'ORF2SGD'")
        result = cursor.fetchall()
        cursor.close()

        if result == []:
            return False
        else:
            return True

    def build_SGDID_map_db(self):
        """
        This adds ORF2SGD table into SGD database
        """
        temp_ORF_set = set()
        cursor = self.SGD_DB_connection.cursor()
        cursor.execute('CREATE TABLE ORF2SGD (ORF VARCHAR(30), SGDID VARCHAR(30))')

        with open('dbxref.tab','r') as f:
            reader = csv.reader(f, delimiter='\t')
            for row in reader:
                if not {row[3]}.issubset(temp_ORF_set):
                    temp_ORF_set.add(row[3])
                    cursor.execute("INSERT INTO ORF2SGD (ORF, SGDID) VALUES(%s, %s)", row[3:5])

        self.SGD_DB_connection.commit()
        cursor.close()

    #TODO: no need to this anymore?!
    def convert_to_SGDID(self, input_ORF):
        cursor = self.SGD_DB_connection.cursor()
        cursor.execute("SELECT SGDID FROM ORF2SGD WHERE ORF = '" + input_ORF + "';")

        ret = cursor.fetchall()
        cursor.close()

        return ret[0][0].encode()


    def set_annot_dict(self): #Payam correct db name
        """
        sets a dictonary that maps annotation to all tables contain that annotation
        """
        #TODO:Maybe input annot_list won't match our naming of Columns. We should consider
        # either giving all annotations to user and ask to select from or publish
        # a Readme that explains possible annotations.

        cursor = self.SGD_DB_connection.cursor()
        ###
        # following line returns a tuple per each table, the first element of
        # each tuple is Table name, others are columns' name in that table
        ###
        query = "SELECT table_name,GROUP_CONCAT(column_name ORDER BY ordinal_position) FROM information_schema.columns WHERE table_schema = '" + self.SGD_db + "' GROUP BY table_name ORDER BY table_name"
        cursor.execute(query)
        tables_columns = cursor.fetchall()
        cursor.close()

        annot2Table = {}
        for t in tables_columns:
            currTable = t[0].encode()
            currAnnots = t[1].encode().split(',')

            for a in currAnnots:
                if a in annot2Table.keys():
                    annot2Table[a].append(currTable)

                else:
                    annot2Table[a] = []
                    annot2Table[a].append(currTable)

        self.annot_dict = annot2Table
        all_annotations = self.annot_dict.keys()

        return all_annotations


    def check_annot(self, annot):
        exists = set(annot).issubset(set(self.annot_dic.keys()))
        return exists

    def do_query(self):
        cursor = self.SGD_DB_connection.cursor()

        all_query_res = []
        for cur_query in self.annlist:
            tables = self.annot_dict[cur_query]
            for curr_table in tables: # Anyway It's wierd to have multiple tables map to the same annotation!
                for curr_gene in self.glist:
                    homologs_SGDID = self.ortho_map[curr_gene]['SGDID']


                    # a SGDID column should be added to all tables in SGD1
                    my_query = " WHERE SGDID IN " +  self.gene_list_to_char(homologs_SGDID) +";"
                    cursor.execute("SELECT " + "SGDID, " + cur_query + " FROM " + curr_table + my_query)
                    # The output format need to be improved
                    all_query_res.append(cursor.fetchall())

        cursor.close()
        return all_query_res
