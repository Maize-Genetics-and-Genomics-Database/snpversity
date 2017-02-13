#!/usr/bin/python
import psycopg2, re,json
#To be used with Tassel-generated "Positions" json file based of Allzeagbsv27 datasets
def parsejson(file,cursor):
    fo = open(file).read()
    cursor.execute("SELECT max(pos) FROM b73v2models")
    max_pos = cursor.fetchall()
    cursor.execute("SELECT max(chr) FROM b73v2models")
    max_chr = cursor.fetchall()
    data = json.loads(fo)
    print("Loaded dataset.")
    transcripts = get01transcripts(cursor)
    print "Loaded transcripts."
    start = False
    for obj in data["PositionList"]:
        chrom = obj["chr"]["name"]
        pos = obj["position"]
        if (chrom >= max_chr[0][0] and pos >= max_pos[0][0]) or start:
            start = True
            tuple_list = getGeneModels(chrom, pos, cursor,transcripts)
            populatedb(chrom,pos,tuple_list,cursor)
        else:
            print "skipped Chr "+str(chrom)+" " + str(pos)

def populatedb(chrom, pos, model_tuple, cursor):
    """
    adds entry to DB
    :param chrom: chromosome
    :param pos: site
    :param cursor: db-connector
    :return: updated value in DB
    """
    if (model_tuple):
        for (category, model) in model_tuple:
            query_populate = "INSERT INTO b73v2models VALUES " \
                         " ('"+category+"', '" + model+"', "+str(chrom)+", " + str(pos) + ");"
            cursor.execute(query_populate)
    else:
        query_populate = "INSERT INTO b73v2models VALUES " \
                         " ('IGR', '', "+str(chrom)+", " + str(pos) + ");"
        cursor.execute(query_populate)


def get01transcripts(cursor):
    """Filters out all 01-transcript genes with corresponding ID"""
    cursor.execute("""SELECT description FROM gene_models WHERE model='mRNA'""")
    rows = cursor.fetchall()
    transcriptions = []
    for row in rows:
        cols = re.split(';', row[0])
        temp = (cols[0][5:],cols[2][3:]) #type/category, transcript-id
        if re.split('_',temp[0])[1] == "T01":
            transcriptions.append(temp)
            #print("added "+temp[0]+" => "+temp[1])
    return transcriptions


def getGeneModels(chrom, pos, cursor,transcripts):
    """
    :return: A list of tuples [intron or exon or IGR, Gene]
    """
    gene_models = []
    if (chrom == 0):
        chr_full = "UNMAPPED"
    else:
        chr_full = "Chr" + str(chrom)

    query_range = "SELECT model,description FROM gene_models " +\
                "WHERE chr='" + chr_full +\
                "' AND " + str(pos) + " >= starts "+\
                "AND " + str(pos) + " < ends AND model != 'mRNA' AND model != 'gene';"
    cursor.execute(query_range)
    matches_range = cursor.fetchall()
    for row in matches_range:
        duplicate = False
        parent_full = re.split(';',row[1])[1] #get Parent=65923
        parent_id = re.split('=',parent_full)[1] #65923
        for (model, id_no) in transcripts:
            if parent_id == id_no:
                for (type, model_prev) in gene_models: #check if gene was encountered previously
                    if model == model_prev:
                        duplicate = True
                    else:
                        duplicate = False
                if (not duplicate):
                    temp = (row[0], model) #type, gene-model
                    print "chr "+str(chrom) + " site " + str(pos) + " with parent-gene " + model + " gives " + row[0]
                    gene_models.append(temp)
    return gene_models


#Main
try:
    conn = psycopg2.connect(database="postgres", user="postgres", password="david101", host="127.0.0.1", port="5432")
    conn.autocommit = True
except:
    print("Could not connect to DB")

cur = conn.cursor()
#getGeneModels(1, 35655673, cur)
#tuple_example = getGeneModels(1,35655673,cur)
#populatedb(1,35655673,tuple_example,cur)
parsejson("v2pos.json",cur)
