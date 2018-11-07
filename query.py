import csv
import mysql.connector as connector
import os
import sys
import subprocess

class Query:
    def __init__(self, species, gene_list, annot_list, dbName, db_user_name, db_password):
        self.spc = species
        self.glist = gene_list
        self.annlist = annot_list
        self.spc_table = species + '_' + 'SC'
        self.db = dbName
        self.user_name = db_user_name
        self.password = db_password
        
    def get_glist(self):
        return(self.glist)
    
    def get_annlist(self):
        return(self.annlist)
    
    def get_spc(self):
        return(self.spc)
    
    def gene_list_to_char(self, list_of_genes):
        my_char="("
        for gene in list_of_genes:
            my_char = my_char + "'" + gene + "'" + ", "
        my_char = my_char[0:(len(my_char)-2)]
        my_char = my_char + ")"
        return(my_char)
    
    
    def encode_to_string(self, encode_tuple_list): # This function can be removed (replaced where it's referenced), or generalized, the current format is too spcific
        # This function converts the tuple list output of fetchall() to a list of tuples of str
        aa_char_list = []
        for x in encode_tuple_list:
            aa_char_list.append((x[0].encode(), x[1].encode()))
        return(aa_char_list)
    
    def set_all_genes(self):
        # this function gets the list of genes corresponding to the spicies in the orhology table
        database_connection = connector.connect(host="localhost",
                                                user=self.user_name,
                                                passwd=self.password,
                                                database=self.db)
        cursor = database_connection.cursor()
        cursor.execute("SELECT " + self.get_spc() + '_genes' + " FROM " + self.spc_table + ";")
        gname_encode = cursor.fetchall()
        database_connection.commit()
        cursor.close()
        gname_str = []
        for gtuple in gname_encode:
            gname_str.append(gtuple[0].encode())
        self.all_genes = gname_str
        
    def get_all_genes(self):
        return(self.all_genes)
    
    def check_gene_exist(self):
        check=set(self.get_glist()).issubset(self.get_all_genes())
        return(check)
    
    def set_orthopairs(self):
        """
        This function retrieves the orthologs of the gene_list in SC
        It reads the orthologs from a table in the database. the table is named: 'desired_species_name'_SC
        and the column is named 'desired_species_name' . Each row represents a pair of orthologs.
        """
        # print "------Connecting to database------"
        database_connection = connector.connect(host="localhost",
                                                user=self.user_name,
                                                passwd=self.password,
                                                database=self.db)
        cursor = database_connection.cursor()
        # print "------connected to " + self.db + "------\n"
        
        ortho_query = " WHERE " + self.get_spc() + '_genes' + " IN " +  self.gene_list_to_char(self.get_glist()) +";"
        cursor.execute("SELECT " + self.get_spc()+ '_genes' + ", SC_genes FROM " + self.spc_table + ortho_query)
        my_ortho_pairs = cursor.fetchall()
        database_connection.commit()
        cursor.close()
        
        if my_ortho_pairs == []:
            print "None of the provided genes exist in the orthology table."
            return(0)
        my_ortho_pairs_str = self.encode_to_string(my_ortho_pairs)
        my_SC_ortho = [x[1] for x in my_ortho_pairs_str]
        self.ortho_pairs_str = my_ortho_pairs_str
    
    def get_orthopairs(self):
        self.set_orthopairs()
        return(self.ortho_pairs_str)
    
    # TODO : detect what is the format of the input gene name : ORF, SGDID or GENE_NAME
    # all gene related tables have to have an SGDID column.
    def set_SGDID(self):
        database_connection = connector.connect(host="localhost",
                                                user=self.user_name,
                                                passwd=self.password,
                                                database=self.db)
        cursor = database_connection.cursor()
        # for now assuming the input is of type ORF
        map_query = "WHERE " + "ORF" + " IN " +  self.gene_list_to_char([x[1] for x in self.get_orthopairs()]) +";"
        cursor.execute("SELECT SGDID FROM SGD_ID_MAP " + map_query)
        SGDid_encode = cursor.fetchall()
        database_connection.commit()
        cursor.close()
        
        SGDid_str = []
        for gtuple in SGDid_encode:
            SGDid_str.append(gtuple[0].encode())
        self.all_genes_SGDid = SGDid_str
    
    def get_SGDID(self):
        return(self.all_genes_SGDid)

    def set_all_annot(self):
        all_possible_annot = []
        database_connection = connector.connect(host="localhost",
                                                user=self.user_name,
                                                passwd=self.password,
                                                database=self.db)
        cursor = database_connection.cursor()
        cursor.execute("SELECT COLUMN_NAME FROM INFORMATION_SCHEMA.COLUMNS WHERE TABLE_SCHEMA='" + self.db + "';")
        all_annot_encode = cursor.fetchall()
        
        database_connection.commit()
        cursor.close()
        for cl in all_annot_encode:
            all_possible_annot.append(cl[0].encode())
        self.all_annot = all_possible_annot
    
    def get_all_annot(self):
        return(self.all_annot)
    
    
    def set_annot_dict(self):
        annot_dict={}
        database_connection = connector.connect(host="localhost",
                                                user=self.user_name,
                                                passwd=self.password,
                                                database=self.db)
        cursor = database_connection.cursor()
        for annot in self.get_all_annot():
            cursor.execute("SELECT DISTINCT TABLE_NAME FROM INFORMATION_SCHEMA.COLUMNS WHERE COLUMN_NAME = '" + annot + "' AND TABLE_SCHEMA='" + self.db + "';")
            my_tables_encode = cursor.fetchall()
            annot_dict[annot] = []
            for tab in my_tables_encode:
                annot_dict[annot].append(tab[0].encode())
        
        database_connection.commit()
        cursor.close()
        self.annot_dic = annot_dict
        
    def get_annot_dict(self):
        return(self.annot_dic)
            
    def check_annot(self):
        check=set(self.get_annlist()).issubset(self.get_all_annot())
        return(check)
    
    def do_query(self):
        database_connection = connector.connect(host="localhost",
                                                user=self.user_name,
                                                passwd=self.password,
                                                database=self.db)
        
        cursor = database_connection.cursor()
        
        if(not self.check_annot()): # check if annotations are present in the database
            print "The provided annotations are not supported in the database."
            return(0)
        
        if(not self.check_gene_exist()): # check if genes are present in the orthology table of the species
            print "The provided genes are not present in the database / Or the naming system is different."
            print "This is an example of the gene names in the orthology table: \n" + self.get_all_genes()[0:5]
            return(0)
        all_query_res = []
        
        for cur_query in self.get_annlist():
            cur_table = self.get_annot_dict()[cur_query][0]
            my_query = " WHERE " + cur_table + ".SGDID" + " IN " +  self.gene_list_to_char(self.get_SGDID()) +";"
            cursor.execute("SELECT " + "SGDID, " + cur_query + " FROM " + cur_table + my_query)
            all_query_res.append(cursor.fetchall())
            
        database_connection.commit()
        cursor.close()
        return(all_query_res)


