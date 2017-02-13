import psycopg2,sys



def getStocksWithInbred(out_file_name, cur, db_table, stock_col, id_col,inbred_col, project=""):
    if project:
        query = "SELECT DISTINCT ON ({0}) {0},{3},{4} FROM {2} WHERE project='{1}';".format(stock_col, project, db_table, inbred_col,id_col)
    else:
        query = "SELECT DISTINCT ON ({0},project) {0},{1},{3} FROM {2};".format(stock_col, inbred_col, db_table,id_col)
    cur.execute(query)
    rows = cur.fetchall()  # (0,1) => (Stock,inbred)
    with open(out_file_name, 'w') as fp:
        for row in rows:
            if row[2]:  # Inbred
                fp.write(row[2] + '\n')
            else:
                fp.write(row[0] + ":" + row[1] + '\n')
        fp.close()
    print "Wrote file {0} containing {1} stocks.".format(out_file_name,str(len(rows)))


def main():
    # Main
    try:
        conn = psycopg2.connect(database="postgres", user="postgres", password="david101", host="127.0.0.1", port="5432")
        conn.autocommit = True
    except:
        print("Could not connect to DB")

    cur = conn.cursor()
    project = ""
    if len(sys.argv) == 2:
        out_file = sys.argv[1]
    elif len(sys.argv) >= 3:
        out_file = sys.argv[1]
        project = ' '.join(sys.argv[2:])
    else:
        sys.exit("Please give output file.")

    db_table = "allzeagbsv27"
    stock_col = "dna_sample"
    inbred_col = "inbred"
    id_col = "lib_prep_id"
    getStocksWithInbred(out_file,cur,db_table,stock_col,inbred_col,id_col,project)

if __name__ == "__main__":
    main()
