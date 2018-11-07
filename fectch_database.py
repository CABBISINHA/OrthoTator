import csv
import requests
import io
import mysql.connector as connector


def get_url_content(url):
    response =requests.get(url).content
    cr = csv.reader(io.StringIO(response.decode('utf-8')), delimiter='\t')
    return cr


def get_dict():
    dict ={}
    with open("dbxref.tab", 'r') as f:
        lines = f.readlines()
        for line in lines:
            line = line.strip('\n')
            temp = line.split('\t')
            key = temp[-1]
            if key is '':
                continue
            val = temp[-2]
            dict[key]=val
    return dict


def connect_to_databse():
    database_connection = connector.connect(host="localhost",user="root",passwd="wjt98")
    cursor = database_connection.cursor()
    cursor.execute("CREATE DATABASE SGD")
    return cursor


def write_to_database_pathway():
    ##TODO: DONE!!!!!!!!!!!!!!!!
    url = "https://downloads.yeastgenome.org/curation/literature/biochemical_pathways.tab"
    csv = get_url_content(url)
    mapping = get_dict()
    data = []
    for n in csv:
        if n is '':
            data.append(n+[""])
        elif n[3].upper() not in mapping:
            data.append(n+[""])
        else:
            data.append(n+[mapping[n[3].upper()]])
    database_connection = connector.connect(host="localhost",user="root",passwd="wjt98")
    cursor = database_connection.cursor()
    cursor.execute("USE SGD")
    cursor.execute("DROP TABLE biochemical_pathways")
    cursor.execute("CREATE TABLE biochemical_pathways (biochem_pathway VARCHAR(100), enzyme VARCHAR(100), EC_number VARCHAR(10), gene_name VARCHAR(10), "
                   "reference VARCHAR(300), SGDID VARCHAR(10))")
    for row in data:
        cursor.execute('INSERT INTO biochemical_pathways(biochem_pathway, enzyme, EC_number, gene_name, reference , SGDID) VALUES(%s, %s, %s, %s, %s, %s)', row[0:6])
    print("Done writing.")
    database_connection.commit()
    cursor.close()


def write_to_database_gene_association():
    ##TODO: DONE!!!!!!!!!!!!!!!!
    database_connection = connector.connect(host="localhost",user="root",passwd="wjt98")
    cursor = database_connection.cursor()
    data = []
    with open ("gene_association.sgd",'r') as f:
        lines = f.readlines()[8:]
        for line in lines:
            data.append(line.strip().split('\t'))
    cursor.execute("USE SGD")
    cursor.execute("DROP TABLE gene_association")
    cursor.execute('CREATE TABLE gene_association (DB VARCHAR(5), SGDID VARCHAR(10), '
                   'Symbol VARCHAR(20), Qualifier VARCHAR(35), GO_Id VARCHAR(12), DB_Reference VARCHAR(100), '
                   'Evidence VARCHAR(5), With_or_From VARCHAR(550), Aspect VARCHAR(5), GeneProductName VARCHAR(100), '
                   'Synonym VARCHAR(200), Type VARCHAR(6), Taxon VARCHAR(13), Date DATE, Annotation_Source VARCHAR(15))')
    for row in data:
         cursor.execute('INSERT INTO gene_association(DB , SGDID, Symbol , Qualifier , GO_Id , DB_Reference , Evidence , With_or_From, Aspect ,'
                        'GeneProductName , Synonym , Type , Taxon , Date , Annotation_Source) '
                        'VALUES(%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)', row[0:15])
    print("Done writing.")
    database_connection.commit()
    cursor.close()


def write_to_database_gene_literature():
    ##TODO: DONE!!!!!!!!!!!!!!!!
    url = "https://downloads.yeastgenome.org/curation/literature/gene_literature.tab"
    data = get_url_content(url)
    database_connection = connector.connect(host="localhost",user="root",passwd="wjt98")
    cursor = database_connection.cursor()
    cursor.execute("USE SGD")
    cursor.execute("DROP TABLE gene_literature")
    cursor.execute('CREATE TABLE gene_literature (PubMed_ID VARCHAR(9), Citation VARCHAR(500), '
                   'Gene_name VARCHAR(20), Feature VARCHAR(35), literature_topic VARCHAR(450), SGDID VARCHAR(10))')
    for row in data:
         cursor.execute("INSERT INTO gene_literature(PubMed_ID , Citation , Gene_name , Feature , Literature_topic, SGDID) VALUES(%s, %s, %s, %s, %s, %s)", row[0:6])
    database_connection.commit()
    database_connection.commit()
    cursor.close()


def write_to_database_go_slim_mapping():
    ##TODO: DONE!!!!!!!!!!!!!!!!
    url = "https://downloads.yeastgenome.org/curation/literature/go_slim_mapping.tab"
    data = get_url_content(url)
    database_connection = connector.connect(host="localhost",user="root",passwd="wjt98")
    cursor = database_connection.cursor()
    cursor.execute("USE SGD")
    cursor.execute('DROP TABLE go_slim_mapping')
    cursor.execute('CREATE TABLE go_slim_mapping (ORF VARCHAR(15), Gene VARCHAR(15), SGDID VARCHAR(10), GO_Aspect VARCHAR(1), GO_Slim_term VARCHAR(100), GOID VARCHAR(10), Feature_type VARCHAR(30))')
    for row in data:
        cursor.execute("INSERT INTO go_slim_mapping(ORF , Gene , SGDID , GO_Aspect , GO_Slim_term, GOID, Feature_type) VALUES(%s, %s, %s, %s, %s, %s, %s)", row[0:7])
    database_connection.commit()
    cursor.close()


def write_to_database_interaction_data():
    ##TODO: DONE!!!!!!!!!!!!!!!!
    url = "https://downloads.yeastgenome.org/curation/literature/interaction_data.tab"
    csv = get_url_content(url)
    data = []
    for n in csv:
        ids = n[10].split('|')
        sgd = ids[0].split(':')[1]
        pmid = ids[1].split(':')[1]
        data.append(n[:10]+[sgd] + [pmid]+n[11:])
    print("Done getting data.")
    database_connection = connector.connect(host="localhost",user="root",passwd="wjt98")
    cursor = database_connection.cursor()
    cursor.execute("USE SGD")
    cursor.execute('DROP TABLE interaction_data')
    cursor.execute("CREATE TABLE interaction_data (Feature_Name_Bait VARCHAR(10), Std_Gene_Name_Bait VARCHAR(10), Feature_name_Hit VARCHAR(10), Std_Gene_Name_Hit VARCHAR(50), "
                   "Experiment_Type VARCHAR(50), Genetic_Physical_Interaction VARCHAR(50), Source VARCHAR(50), Manually_curated VARCHAR(50), Notes VARCHAR(400),"
                   "Phenotype VARCHAR(50), SGDID VARCHAR(10), PubMed_ID VARCHAR(10), Citation VARCHAR(500))")
    for row in data:
        cursor.execute("INSERT INTO interaction_data (Feature_Name_Bait, Std_Gene_Name_Bait, Feature_name_Hit, Std_Gene_Name_Hit, Experiment_Type, "
                       "Genetic_Physical_Interaction, Source, Manually_curated, Notes, Phenotype, SGDID, PubMed_ID, Citation) VALUES"
                       "(%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)",row[0:13])
    print("Done writing.")
    database_connection.commit()
    cursor.close()


def write_to_database_phenotype_data():
    ##TODO: DONE!!!!!!!!!!!!!!!!
    url = "https://downloads.yeastgenome.org/curation/literature/phenotype_data.tab"
    data = get_url_content(url)
    database_connection = connector.connect(host="localhost",user="root",passwd="wjt98")
    cursor = database_connection.cursor()
    cursor.execute("USE SGD")
    cursor.execute('DROP TABLE phenotype_data')
    cursor.execute('CREATE TABLE phenotype_data (Feature_Name VARCHAR(10), Feature_Type VARCHAR(50), Gene_Name VARCHAR(10), SGDID VARCHAR(10), '
                   'Reference VARCHAR(40), Experiment_Type VARCHAR(300), Mutant_Type VARCHAR(50), Allele VARCHAR(200), Strain_Background VARCHAR(100), '
                   'Phenotype VARCHAR(100), Chemical VARCHAR(200), Observed_condition VARCHAR(200), Details VARCHAR(400), Reporter VARCHAR(200))')

    for row in data:
        #print(row)
        cursor.execute("INSERT INTO phenotype_data (Feature_Name, Feature_Type, Gene_Name, SGDID, Reference, Experiment_Type,"
                       "Mutant_Type, Allele, Strain_Background, Phenotype, Chemical, Observed_condition, Details, Reporter) "
                       "VALUES(%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)",row[0:14])
    database_connection.commit()
    cursor.close()


def write_to_database_sgd_features():
    ##TODO: DONE!!!!!!!!!!!!!!!!
    url = "https://downloads.yeastgenome.org/curation/chromosomal_feature/SGD_features.tab"
    data = get_url_content(url)
    database_connection = connector.connect(host="localhost",user="root",passwd="wjt98")
    cursor = database_connection.cursor()
    cursor.execute("USE SGD")
    cursor.execute('DROP TABLE SGD_features')
    cursor.execute('CREATE TABLE SGD_features (SGDID VARCHAR(10), Feature_Type VARCHAR(50), Feature_qualifier VARCHAR(50), Feature_name VARCHAR(50), '
                   'Standard_gene_name VARCHAR(50), Alias VARCHAR(200), Parent_feature_name VARCHAR(20), 2nd_SGDID VARCHAR(100), Chromosome VARCHAR(30), '
                   'Start_coordinate VARCHAR(20), Stop_coordinate VARCHAR(20), Strand VARCHAR(20), Genetic_position VARCHAR(50), '
                   'Coordinate_version VARCHAR(50), Sequence_version VARCHAR(50), Description VARCHAR(500))')
    for row in data:
        cursor.execute('INSERT INTO SGD_features (SGDID, Feature_Type, Feature_qualifier, Feature_name, Standard_gene_name, '
                       'Alias, Parent_feature_name, 2nd_SGDID, Chromosome, Start_coordinate, Stop_coordinate, Strand, Genetic_position, '
                       'Coordinate_version, Sequence_version, Description) VALUES(%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)',row[0:17])
    database_connection.commit()
    cursor.close()


if __name__ == "__main__":
    #connect_to_databse()
    #write_to_database_pathway()
    #write_to_database_gene_association()
    #write_to_database_gene_literature()
    #write_to_database_go_slim_mapping()
    write_to_database_interaction_data()
    #write_to_database_phenotype_data()
    #write_to_database_sgd_features()
