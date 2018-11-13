import csv
import mysql.connector as connector
import os
import sys
import subprocess
import pandas as pd

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
        if not self.table_exists('ORF2SGD'):
            self.build_SGDID_map_db()
            print "done_building_SGDID_Map"
        self.set_all_genes()
        print "done_setting_all_genes"
        self.all_annots = self.set_annot_dict() # This also sets annot_dict
        print "done_setting_all_annots"
        self.set_ortho_map()
        print "done_setting_ortho_map"
        # Check if all requested annotations are in the database
        for c_annot in self.annlist:
            if(not c_annot in self.all_annots):
                print " annotation " + c_annot +  " not present in the database "
                print " available annotations: "
                print set(self.all_annots)
                #return
                sys.exit()
        # Check if all the requested genes are in the database, removing genes that are not present
        f_gene_list = list(self.glist)
        for c_gene in self.glist:
            if(not c_gene in self.all_genes):
                print "Gene named " +c_gene + " is not present in the orthology database. This gene will be removed from the query."
                f_gene_list.remove(c_gene)

        if(len(f_gene_list) == 0):
            print ("none of the query genes were present in orthology database.")
            print " available genes: "
            print self.all_genes
            sys.exit()
        self.glist = f_gene_list
        #print self.glist


    def gene_list_to_char(self, list_of_genes):
        my_char="("
        for gene in list_of_genes:
            my_char = my_char + "'" + gene + "'" + ", "
        my_char = my_char[0:-2]
        my_char = my_char + ")"
        return my_char


    def set_all_genes(self):
        # this function sets the list of genes corresponding to the species in the orhology table
        cursor = self.inparanoid_DB_connection.cursor()
        cursor.execute("SELECT " + self.spc + '_genes' + " FROM " + self.spc_table + ";")
        gname_encode = cursor.fetchall()
        cursor.close()
        gname_str = []
        for gtuple in gname_encode:
            gname_str.append(gtuple[0].encode())
        self.all_genes = set(gname_str)

    def set_ortho_map(self):
        ## consider adding a threshold as input to filter orthologs if needed
        """
        This function retrieves the orthologs of the gene_list in SC.
        It reads the orthologs from a table in the database. the table is named: 'desired_species_name'_SC
        and the column is named 'desired_species_name' + '_genes'.

        It forms a dictonary, keys are genes in gene_list (those that have orthologs), value is another
        dictopnary with two keys (for now) 'ORF' and 'SGDID', value for each is an array containing all
        orthologs.
        """
        cursor = self.inparanoid_DB_connection.cursor()

        ortho_query = " WHERE " + self.spc + '_genes' + " IN " +  self.gene_list_to_char(self.glist) +";"
        cursor.execute("SELECT " + self.spc+ '_genes' + ", SC_genes FROM " + self.spc_table + ortho_query)
        my_ortho_pairs = cursor.fetchall()
        #print my_ortho_pairs
        cursor.close()
        if my_ortho_pairs == []:
            print "None of the provided genes exist in the orthology table."
            return
        orthoMap = {}
        my_all_SGD = []
        for p in my_ortho_pairs:
            key = p[0].encode()
            my_all_SGD.append(self.convert_to_SGDID(p[1].encode()))
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
        self.SGD_set = set(my_all_SGD)


    def table_exists(self, table_name):
        cursor = self.SGD_DB_connection.cursor()
        cursor.execute("SHOW TABLES LIKE " + "'" + table_name + "'")
        result = cursor.fetchall()
        cursor.close()
        if result == []:
            return False
        else:
            return True


    def table_drop(self, table_name):
        my_cursor = self.SGD_DB_connection.cursor()
        my_cursor.execute("SHOW TABLES LIKE " + "'" + table_name + "'")
        result = my_cursor.fetchall()
        if not result == []:
            my_cursor.execute("DROP TABLE " + table_name )
        my_cursor.close()


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

    def convert_to_ORF(self, input_SGDID):
        cursor = self.SGD_DB_connection.cursor()
        cursor.execute("SELECT ORF FROM ORF2SGD WHERE SGDID = '" + input_SGDID + "';")
        ret = cursor.fetchall()
        cursor.close()
        return ret[0][0].encode()

    def create_cur_qu2ORF_table(self):
        """
        creates a temporary table 'cur_qu2ORF_table' that contains the query genes and their corresponding ortholog ORF
        """
        self.table_drop('cur_qu2ORF_table')
        ortho_query = " WHERE " + self.spc + '_genes' + " IN " +  self.gene_list_to_char(self.glist) +";"
        new_cursor = self.SGD_DB_connection.cursor()
        new_cursor.execute("CREATE TABLE cur_qu2ORF_table SELECT " + self.spc + '_genes' + ", SC_genes FROM " + self.inparanoid_db + "." + self.spc_table + ortho_query)
        print "cur_qu2ORF_table created"
        new_cursor.close()

    def create_qgene_ORF_SGDID_Description_table(self):
        """
        Creates a temporary table where the columns are:  query genes, ortholog ORFs, ortholog SGDID and ortholog Description
        """
        self.table_drop('cur_qu2ORF2SGD_table')
        self.table_drop('cur_qu2ORF2SGD_description_table')
        cursor = self.SGD_DB_connection.cursor()
        cursor.execute("CREATE TABLE cur_qu2ORF2SGD_table SELECT " + (self.spc+ '_genes') + ', ORF, '+ 'SGDID' + " FROM  cur_qu2ORF_table, ORF2SGD  WHERE cur_qu2ORF_table.SC_genes = ORF2SGD.ORF")
        cursor.execute("CREATE TABLE cur_qu2ORF2SGD_description_table SELECT " + (self.spc+ '_genes') + ', ORF, '+ 'cur_qu2ORF2SGD_table.SGDID, ' + "Description" +" FROM  cur_qu2ORF2SGD_table, SGD_features  WHERE cur_qu2ORF2SGD_table.SGDID = SGD_features.SGDID")
        cursor.close()


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


    def do_query(self):
        self.create_cur_qu2ORF_table()
        self.create_qgene_ORF_SGDID_Description_table()
        print "done_create_qgene_ORF_SGDID_table"
        cursor = self.SGD_DB_connection.cursor()
        all_query_res = []
        all_query_conc = []
        for cur_query in self.annlist:
            print "cur_query"
            print cur_query
            print "#######"
            tables = self.annot_dict[cur_query]
            for curr_table in tables: # Anyway It's wierd to have multiple tables map to the same annotation!
                print "curr_table"
                print curr_table
                print "#######"
                self.table_drop('cur_table_query_table')
                my_query = " WHERE SGDID IN " +  self.gene_list_to_char(list(self.SGD_set)) +";"
                # create a temporary table containing the SGDID and cur_query
                cursor.execute("CREATE TABLE cur_table_query_table SELECT " + "SGDID, " + cur_query + " FROM " + curr_table + my_query)
                # select to get name, ORF, SGDID, Description and cur_query
                cursor.execute("SELECT " + self.spc+ '_genes' + ', ORF, ' + 'cur_qu2ORF2SGD_description_table.SGDID, ' + ' Description, ' + cur_query + " FROM " + "cur_qu2ORF2SGD_description_table, cur_table_query_table WHERE cur_qu2ORF2SGD_description_table.SGDID = cur_table_query_table.SGDID;")
                all_query_res.append(cursor.fetchall())
                # create a table where all annotations for the same pair are concata
                self.table_drop('temp_table')
                cursor.execute("CREATE TABLE temp_table " + "SELECT " + self.spc+ '_genes' + ', ORF, ' + 'cur_qu2ORF2SGD_description_table.SGDID, ' + ' Description, ' + cur_query + " FROM " + "cur_qu2ORF2SGD_description_table, cur_table_query_table WHERE cur_qu2ORF2SGD_description_table.SGDID = cur_table_query_table.SGDID;")
                cursor.execute("SELECT " + self.spc+ '_genes' + ', ORF, ' + 'SGDID, ' + ' Description, ' + "GROUP_CONCAT("+ cur_query + ")" + " FROM temp_table GROUP BY " + self.spc+ '_genes, ' + ' ORF, ' + 'SGDID, ' + ' Description ;')
                all_query_conc.append(cursor.fetchall())
                self.table_drop('temp_table')


        cursor.execute("DROP TABLE cur_table_query_table")
        cursor.execute("DROP TABLE cur_qu2ORF2SGD_table")
        cursor.execute("DROP TABLE cur_qu2ORF_table")
        cursor.execute("DROP TABLE cur_qu2ORF2SGD_description_table")
        cursor.close()
        return all_query_res, all_query_conc

    def do_query_print(self):
        """
        This is an alternative method to do_query, it's particulary useful for presentation purposes.
        Only uses SGD-DB tables and ortholog tables. Consider merging/replacing with do_query.

        Returns a pandas dataframe.

        In cases of having annotations appearing in multiple tables, for now it just fetchs one of them.
        This may raise error in sql syntax, this should be fixed by considering all tables, and also adding
        a GROUP_BY command to keep track of tables --> This is particulary important if we decided to replace this
        with do_query() function

        Note that it summarizes description for presentation purposes. Make sure to fix it
        if we decided to replace with do_query()
        """
        cursor = self.SGD_DB_connection.cursor()

        join_statement = "INNER JOIN "
        annots_to_fetch_char = ""
        tablesSet = set()
        columns = ['SGDID', 'Description']
        for query in self.annlist:
            if not {self.annot_dict[query][0]}.issubset(tablesSet) and self.annot_dict[query][0] != "SGD_features":
                tablesSet.add(self.annot_dict[query][0]) # for now just 1 table is considered
                join_statement += self.annot_dict[query][0]
                join_statement += " ON " + self.annot_dict[query][0] + ".SGDID = "
                join_statement += "SGD_features.SGDID INNER JOIN "

            annots_to_fetch_char += query
            annots_to_fetch_char += ", "
            columns.append(query)
        join_statement = join_statement[:-12]
        annots_to_fetch_char = annots_to_fetch_char[:-2]
        ret = [] # return a list of pandas data frame, 1 per each gene in gene_list

        columns.insert(0,'Query Gene')
        output_as_array = []
        for curr_gene in self.glist:
            homologs_SGDID = self.ortho_map[curr_gene]['SGDID']
            my_query = " WHERE SGD_features.SGDID IN " +  self.gene_list_to_char(homologs_SGDID) +";"
            cursor.execute('SELECT SGD_features.SGDID, Description, ' + annots_to_fetch_char + " FROM SGD_features " + join_statement + my_query)
            output = cursor.fetchall()

            for line in output:
                currLine = []
                for t in line:
                    if len(t.encode()) > 20:
                        t = t[0:15] + " ..."
                    currLine.append(t)
                currLine[0] = self.convert_to_ORF(currLine[0])
                #currLine[1] = currLine[1][0:20] + " ..."
                currLine.insert(0, curr_gene[0:15]+" ...")
                output_as_array.append(currLine)
        columns[1] = 'ORF'
        df = pd.DataFrame(data = output_as_array, columns = columns)
        cursor.close()
        return df
