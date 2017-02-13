#!/usr/bin/python
import psycopg2, re,json,sys
"""To be used with Tassel-generated "Positions" json file based of HapMapV3"""


def parsejson(file,cursor):
    ranges = get01transcriptranges(cursor)
    print "Loaded transcripts."
    cursor.execute("SELECT max(pos) FROM b73v3models")
    max_pos = cursor.fetchall()
    cursor.execute("SELECT max(chr) FROM b73v3models")
    max_chr = cursor.fetchall()
    chrom = ""
    pos = 0
    chrom_now = 0
    fo = open(file)
    while (True): # Check starting at 1
        ch = fo.read(1)
        if (ch ==''): break
        if (ch == '['):
            chunk = ""
            special_char_count = 0
            special_char = ','
            i = 1
            first_iteration = True
            while (True):#begin json-obj
                ch = fo.read(1)
                if (ch == special_char): #{chr:...
                    special_char_count += 1
                    if (special_char_count == 5): #found json object
                        special_char_count = 0
                        if (chunk[0] ==special_char):
                            chunk = chunk[1:] #if begins with , drop the ,
                        json_string = chunk #i = 234 on first try
                        obj = json.loads(json_string)
                        chrom = int(obj["chr"]["name"])
                        pos = int(obj["position"]) #begin populating DB
                        json_string = ""
                        chunk = ""
                        if not(chrom == chrom_now):
                            chrom_now = chrom
                            exons = getAllExon(chrom,cur)
                        if (chrom >= max_chr[0][0] and pos >= max_pos[0][0]):
                            #print "begin " + str(chrom)+", "+ str(pos)
                            tuple_list = getGeneModels(chrom, pos, cursor,ranges,exons)
                            populatedb(chrom,pos,tuple_list,cursor)
                        else:
                            print "skipped Chr "+str(chrom)+" " + str(pos)
                        i = i+1
                chunk = chunk + ch

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
            query_populate = "INSERT INTO b73v3models VALUES " \
                         " ('"+category+"', '" + model+"', "+str(chrom)+", " + str(pos) + ");"
            cursor.execute(query_populate)
    else:
        query_populate = "INSERT INTO b73v3models VALUES " \
                         " ('IGR', '', "+str(chrom)+", " + str(pos) + ");"
        cursor.execute(query_populate)

def get01transcriptranges(cursor):
    """Filters out all 01-transcript genes with corresponding range
    returns ranges[gene-ID] => (name,start,end,chromosome)"""
    cursor.execute("""SELECT description,starts,ends,chr FROM gene_modelsv3 WHERE model='mRNA'
    OR model='gene' OR model='miRNA_gene' OR model='miRNA'""")
    rows = cursor.fetchall()
    ranges = {}
    for row in rows:
        descriptors = re.split(';', row[0])#[ID=transcript:GRGRMASD312_T01;parent=,...]
        temp = (descriptors[0][5:],descriptors[2][3:]) #drop name=, ID=
       # temp = (int(row[1]),int(row[2]), int(row[3])) #(10045,10051,1) = start,pos,chr
        try:
            ending = re.split('_',temp[0])[1]
            if ending == "T01" or ending == "P01":
                print "added "+temp[1]+" => ("+temp[0]+", "+str(row[1])+", "+str(row[2])+", "+row[3]+")"
                ranges[temp[1]] = (temp[0],row[1],row[2],row[3]) #ranges[gene-ID] => (name,start,end,chromosome)
        except:
            continue #do nothing. (not T01 or P01)
    return ranges



def checkGeneExistence(transcription, gene_models):
    for (type, model_prev) in gene_models: #check if gene was encountered previously
        if transcription == model_prev:
            return True
    return False

def checkIntron(chrom, pos, ranges):
    gene_models = []
    for (transcription_id, temp) in ranges.iteritems():
        duplicate = checkGeneExistence(temp[0], gene_models)
        if (chrom == temp[3] and pos >= temp[1] and pos < temp[2] and not duplicate): #temp = (start,end,chromosome)
            tuple = ("intron", temp[0])
            gene_models.append(tuple) #type, gene-transcript
    return gene_models


def getAllExon(chrom,cursor):
    chrom_full = "Chr" + str(chrom)
    query_range = "SELECT model,description,starts,ends,chr FROM gene_modelsv3 " +\
                "WHERE chr='" + chrom_full+"' AND model != 'miRNA_gene' AND model != 'mRNA' AND model != 'gene' AND model != 'miRNA';"
    cursor.execute(query_range)
    return cursor.fetchall()



def checkExon(chrom, pos, cursor,chrom_exons,ranges):
    gene_models = []
    for row in chrom_exons:
        if (row[4] == chrom and row[2] <= pos and row[3] > pos): #matches
            duplicate = False
            try:
                name_full = re.split(';',row[1])[0] #get name=GRMZM2343_P01
                name = name_full[5:] #name=GRMZM2343_P01
                ending = re.split('_',name)[1] #P01
                if (ending == "P01" or ending=="T01"):
                    duplicate = checkGeneExistence(name, gene_models)
                    if not duplicate:
                        gene_type = row[0]
                        if (gene_type == "CDS"):
                            gene_type = 'exon'
                        tuple = (gene_type, name)
                        gene_models.append(tuple) #type, gene-transcript
            except:#could be T01
                parent_full = re.split(';',row[1])[0] #get Parent=874531
                parent_key = re.split('=',parent_full)[1] #87453
                for (transcript_id, range_tuple) in ranges.iteritems():
                    if parent_key == transcript_id: #for v3 must be exon
                        duplicate = checkGeneExistence(range_tuple[0], gene_models)
                        if not duplicate:
                            gene_type = row[0]
                            if (gene_type == "CDS"):
                                gene_type = 'exon'
                            tuple = (gene_type, range_tuple[0])
                            gene_models.append(tuple) #type, gene-transcript
    return gene_models

def getGeneModels(chrom, pos, cursor, ranges,exons):
    """
    :return: A list of tuples [intron or exon or IGR, Gene]
    """
    #exon = 10299. Intron = 11538
    chrom_full = "Chr" + str(chrom)

    gene_models = checkExon(chrom_full,pos,cursor,exons,ranges)
    if not gene_models:
        gene_models = checkIntron(chrom_full,pos,ranges)
    return gene_models

#Main
try:
    conn = psycopg2.connect(database="postgres", user="postgres", password="david101", host="127.0.0.1", port="5432")
    conn.autocommit = True
except:
    print("Could not connect to DB")

cur = conn.cursor()
#ranges = get01transcriptranges(cur)
parsejson("v3pos.json",cur)
