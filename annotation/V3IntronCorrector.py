import psycopg2


def correctIntrons(cur, chrom):
    step = 10000
    start = 0
    end = 301410279
    for offset in range(start,end,step):
        models = getModels(cur,step,offset,chrom)
        ranges = getFalseIntronPos(cur,models)
        updateDB(ranges,cur)
        end = offset
    print "FINISHED correcting chromosome "+str(chrom)


def getFalseIntronPos(cur,models):
    ranges = []
    for i in range(len(models)): #i = start
        try:
            if models[i][2] == "IGR" and models[i+1][2] == "intron" and models[i+2][2] == "IGR":
                temp = (models[i+1][1], models[i+1][3], models[i+1][0])#start,end,chr
                ranges.append(temp)
        except:
            print "EXCEPTION @ chr " + str(models[i][0])+" pos "+str(models[i][1])
    #print "loaded false intron ranges."
    return ranges


def updateDB(ranges,cursor):
    for range in ranges:
        query = "UPDATE b73v3ranges SET type='exon' WHERE chr="+str(range[2])+" AND pos >= "+str(range[0])+" AND pos <= "+str(range[1])+";"
        print query
        #cursor.execute(query)


def getModels(cursor,limit,offset,chrom):
    query = "select chr,pos,type,ends from b73v3ranges WHERE chr="+str(chrom)+" ORDER BY chr,pos LIMIT "+str(limit)+" OFFSET "+str(offset)+";"
    cursor.execute(query)
    #print "loaded models."
    return cursor.fetchall()

# Main
try:
    conn = psycopg2.connect(database="postgres", user="postgres", password="david101", host="127.0.0.1", port="5432")
    conn.autocommit = True
except:
    print("Could not connect to DB")

cur = conn.cursor()
for i in range(1,11):
    correctIntrons(cur, i)