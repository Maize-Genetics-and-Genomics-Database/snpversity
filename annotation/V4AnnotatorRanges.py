#!/usr/bin/python
import psycopg2, re, sys
# Takes in a chromosome index to look through the b73v4models table and generate corresponding b73v4ranges annotations. 
# IMPORTANT: Be sure to run for each chromosome and correct afterwards with V4IntronCorrector.



def annotate(cursor,current_chrom, current_pos, end):
    ranges = get01transcriptranges(cursor)
    next_batch = get_remaining_batch_from_chrom(cursor, current_chrom, current_pos)
    interfeatures(next_batch)

    for r in next_batch:
        if r[1] > end:
            break
        else:
            annotations = get_annotations_from_gene_old(r, cursor,ranges) # [(Type, Gene-name, start, end, chrom), ...]
            if annotations:
                interfeature_introns(annotations)
                populatedbFromAnnotations(annotations,cursor)
    """ models = getAllAnnotationsFromChrom(r[3], cursor)
     tuple_list = getGeneModels(r[3],r[1],cursor,transcripts,models)
     for t in tuple_list:
         print t"""


def get_remaining_batch_from_chrom(cursor, current_chrom, current_pos):
    """Gets all gene positions of T001 or P001 transcripts"""
    query_next_pos = "SELECT * FROM gene_modelsv4 WHERE chr='Chr{0}' AND " \
                     "starts>={1} AND " \
                     "(model='mRNA'OR model='gene' OR model='miRNA_gene' OR model='miRNA')" \
                     " ORDER BY starts, ends".format(str(current_chrom),str(current_pos))
    transcript_positions = []
    cursor.execute(query_next_pos)
    transcripts = cursor.fetchall()
    for row in transcripts:
        descriptors = re.split(';', row[8])  # [ID=transcript:GRGRMASD312_T001;parent=,...]
        temp = (descriptors[0][5:], descriptors[2][3:])  # drop name=, ID=
        # temp = (int(row[1]),int(row[2]), int(row[3])) #(10045,10051,1) = start,pos,chr
        try:
            ending = re.split('_', temp[0])[1]
            if ending == "T001" or ending == "P001":
                ranges = (temp[0], row[3], row[4], current_chrom)  # (name,start,end,chromosome)
                transcript_positions.append(ranges)
        except:
            continue  # do nothing. (not T001 or P001)
    return transcript_positions


def interfeatures(ranges, type=""):
    """Adds IGR's"""
    for i in xrange(0,len(ranges)*2-2,2):
        start_IGR = ranges[i][2]
        end_IGR = ranges[i+1][1]
        chrom_IGR = ranges[i+1][3]
        tuple_IGR = (type, start_IGR, end_IGR, chrom_IGR)
        ranges.insert(i+1, tuple_IGR)

def interfeature_introns(annotations, type_feat='intron'):
    """
    Accepts annotations of type #[(Type, Gene-name, start, end, chrom), ...] Inserts introns in appropriate regions.
    :param annotations:
    :return:
    """
    for i in xrange(0,len(annotations)*2-2,2):
        try:
            if annotations[i][0] in ['exon', 'three_prime_UTR', 'five_prime_UTR'] \
            and annotations[i+1][0] in ['exon', 'three_prime_UTR', 'five_prime_UTR']\
            and annotations[i][1] == annotations[i+1][1]:
                start_intron = annotations[i][3]
                end_intron = annotations[i+1][2]
                name_intron = annotations[i][1]
                chrom_intron = annotations[i+1][4]
                tuple = (type_feat, name_intron, start_intron, end_intron, chrom_intron)
                annotations.insert(i+1, tuple)
        except IndexError:
            print "Error addinng intron to {0}".format(annotations[0][1])
            break


def populatedb(chrom, pos,end, model_tuple, cursor):
    """
    adds entry to DB
    :param chrom: chromosome
    :param pos: site
    :param cursor: db-connector
    :return: updated value in DB
    """
    for (type, model) in model_tuple:
        query_populate = "INSERT INTO b73v4ranges VALUES " \
                     " ('" + type + "', '" + model + "', " + str(chrom) + ", " + str(pos) + ", "+str(end)+");"
        cursor.execute(query_populate)
        print(query_populate)


def populatedbFromAnnotations(annotations, cursor):
    for a in annotations:
        annotation_type, gene_name, start, end, chrom = a
        query_populate = "INSERT INTO b73v4ranges VALUES " \
                     " ('" + annotation_type + "', '" + gene_name + "', " + str(chrom) + ", " + str(start) + ", "+str(end)+");"
        cursor.execute(query_populate)
        print(query_populate)


def get01transcriptranges(cursor):
    """Filters out all 01-transcript genes with corresponding range
    returns ranges[gene-ID] => (name,start,end,chromosome)"""
    cursor.execute("""SELECT description,starts,ends,chr FROM gene_modelsv4 WHERE model='mRNA'
    OR model='gene' OR model='miRNA_gene' OR model='miRNA'""")
    rows = cursor.fetchall()
    ranges = {}
    for row in rows:
        descriptors = re.split(';', row[0]) #[ID=transcript:GRGRMASD312_T001;parent=,...]
        temp = (descriptors[0][5:],descriptors[2][3:]) #drop name=, ID=
       # temp = (int(row[1]),int(row[2]), int(row[3])) #(10045,10051,1) = start,pos,chr
        try:
            ending = re.split('_',temp[0])[1]
            if ending == "T001" or ending == "P001":
                #print "Adding to transcripts:"+temp[1]+" => ("+temp[0]+", "+str(row[1])+", "+str(row[2])+", "+row[3]+")"
                ranges[temp[1]] = (temp[0],row[1],row[2],row[3]) #ranges[gene-ID] => (name,start,end,chromosome)
        except:
            continue #do nothing. (not T001 or P001)
    return ranges


def get_annotations_from_gene(gene_tuple, cur):
    """
    :param gene_tuple:
    :param cur:
    :return: List of annotations of the form (Type, Gene-name, start, end, chrom)
    """
    # get results from gene or IGR
    name, start_gene, end_gene, chrom = gene_tuple
    # IGR
    if not name:
        return [("IGR","",start_gene,end_gene,chrom)]

    query_get_annotations = "SELECT model,description,starts,ends,chr from gene_modelsv4 WHERE chr='Chr{0}' AND " \
                            "starts>={1} AND ends<={2} AND " \
                            "(model='CDS' OR model='three_prime_UTR' OR model='five_prime_UTR')".format(str(chrom),str(start_gene),str(end_gene)) # AND (model='CDS' OR model='three_prime_UTR' OR model='five_prime_UTR'
    cur.execute(query_get_annotations)
    annotations = cur.fetchall()
    # Full Intron
    if not annotations:
        return [("Intron", name, start_gene, end_gene, chrom)]
    gene_annotations_full = []
    prev_end = start_gene
     # Starts with intron
    if annotations[0][2] >= start_gene:
        tuple_next = ("Intron", name, start_gene, annotations[0][2], chrom)
        prev_end = annotations[0][2]
        gene_annotations_full.append(tuple_next)
    for a in annotations:
        gene_type = a[0]
        # Intron between exons
        if a[2] <= prev_end:
            tuple_next = ("Intron", name, prev_end, a[2], chrom)
            gene_annotations_full.append(tuple_next)
        if a[0] == "CDS":
            gene_type = "Exon"
        tuple_next = (gene_type, name, a[2], a[3], chrom)
        gene_annotations_full.append(tuple_next)
        prev_end = a[3]
    return gene_annotations_full

def get_annotations_from_gene_old(gene_tuple, cur, ranges):
    """
    :param gene_tuple:
    :param cur:
    :return: List of annotations of the form (Type, Gene-name, start, end, chrom)
    """
    # get results from gene or IGR
    name, start_gene, end_gene, chrom = gene_tuple
    # IGR
    if not name:
        return [("IGR","",start_gene,end_gene,chrom)]

    query_get_annotations = "SELECT model,description,starts,ends,chr from gene_modelsv4 WHERE chr='Chr{0}' AND " \
                            "starts>={1} AND ends<={2} AND " \
                            "(model='CDS' OR model='three_prime_UTR' OR model='five_prime_UTR') ORDER BY starts, ends".format(str(chrom),str(start_gene),str(end_gene)) # AND (model='CDS' OR model='three_prime_UTR' OR model='five_prime_UTR'
    cur.execute(query_get_annotations)
    annotations = cur.fetchall()
    # Full Intron
    if annotations:
        gene_annotations_full = []
        annotation_start = start_gene
        annotation_end = annotation_start
        for a in annotations:
            gene_models = getGeneModels(chrom,a[2],cur,ranges,annotations)
            for g in gene_models:
                """if a[2] > annotation_end:
                    tuple_next = ("intron", name, annotation_end, a[2], chrom)
                    gene_annotations_full.append(tuple_next)"""
                tuple_next = (g[0], g[1], a[2], a[3], chrom)
                gene_annotations_full.append(tuple_next)
                annotation_end = a[3]
        return gene_annotations_full
    else:
        return [("exon", name, start_gene, end_gene, chrom)]



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


def getAllAnnotationsFromChrom(chrom, cursor):
    chrom_full = "Chr" + str(chrom)
    query_range = "SELECT model,description,starts,ends,chr FROM gene_modelsv4 " +\
                "WHERE chr='" + chrom_full+"' AND model != 'miRNA_gene' AND model != 'mRNA' AND model != 'gene' AND model != 'miRNA';"
    cursor.execute(query_range)
    return cursor.fetchall()



def checkExon(chrom, pos, cursor, annotations, ranges):
    gene_models = []
    for row in annotations:
        if (row[4] == chrom and row[2] <= pos and row[3] > pos): #matches
            duplicate = False
            try:
                name_full = re.split(';',row[1])[0] #get name=GRMZM2343_P001
                name = name_full[5:] #name=GRMZM2343_P001
                ending = re.split('_',name)[1] #P001
                if (ending == "P001" or ending=="T001"):
                    duplicate = checkGeneExistence(name, gene_models)
                    if not duplicate:
                        gene_type = row[0]
                        if (gene_type == "CDS"):
                            gene_type = 'exon'
                        tuple = (gene_type, name)
                        gene_models.append(tuple) #type, gene-transcript
            except:#could be T001
                parent_full = re.split(';',row[1])[0] #get Parent=874531
                parent_key = re.split('=',parent_full)[1] #87453
                for (transcript_id, range_tuple) in ranges.iteritems():
                    if parent_key == transcript_id: #for v4 must be exon
                        duplicate = checkGeneExistence(range_tuple[0], gene_models)
                        if not duplicate:
                            gene_type = row[0]
                            if (gene_type == "CDS"):
                                gene_type = 'exon'
                            tuple = (gene_type, range_tuple[0])
                            gene_models.append(tuple) #type, gene-transcript
    return gene_models

def getGeneModels(chrom, pos, cursor, ranges, annotations):
    """
    :return: A list of tuples [(intron or exon or IGR, Gene-model1),(intron or exon or IGR, Gene-model2),...]
    """
    chrom_full = "Chr" + str(chrom)

    gene_models = checkExon(chrom_full, pos, cursor, annotations, ranges)
    if not gene_models:
        gene_models = checkIntron(chrom_full,pos,ranges)
    if not gene_models:
        gene_models = [("IGR","")]
    return gene_models

#Main
try:
    conn = psycopg2.connect(database="postgres", user="postgres", password="david101", host="127.0.0.1", port="5432")
    conn.autocommit = True
except:
    print("Could not connect to DB")

cur = conn.cursor()
try:
    max_chr = sys.argv[1]
except IndexError:
    sys.exit()
cur.execute("SELECT max(pos) FROM b73v4ranges WHERE chr="+str(max_chr))
max_pos = cur.fetchall()
cur.execute("SELECT max(ends) FROM gene_modelsv4 WHERE chr='Chr{0}'".format(str(max_chr)))
end = cur.fetchall()
#max_pos = [[16167234]]
#end=[[16221585]]
annotate(cur,max_chr,max_pos[0][0], end[0][0])
