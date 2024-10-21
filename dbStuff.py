import psycopg2
from utilities import powerset
import random

from psycopg2 import sql


def getSizeOf(conn, table):
    query = "select count(*) from \"" + table + "\";"
    res=execute_query(conn, query)
    return res[0][0]

def dropIndex(conn, table, sel):
    indexname = table + '_' + sel
    query = "drop index if exists \"" + indexname +"\";"
    execute_query(conn, query)

def dropAllIndex(conn, table):
    query="select indexname from pg_catalog.pg_indexes where tablename = \'" + table + "\';"
    res=execute_query(conn, query)
    for r in res:
        if not r[0].startswith('Key') and not r[0].endswith('_pkey'):
            drop="drop index \"" + r[0] + "\";"
            execute_query(conn, drop)

def generateHashIndex(conn, table, sel):
    indexname=table+'_'+sel
    query = "create index if not exists \"" + indexname + "\" on \"" + table + "\" using hash(" + sel + ");"
    execute_query(conn,query)

def generateIndex(conn, table, sel):
        indexname = table + '_' + sel
        query = "create index if not exists \"" + indexname + "\" on \"" + table + "\"(" + sel + ");"
        execute_query(conn, query)

# generate a multicolumn index on table
def generateMulticolIndex(conn, table, list, selAtt):
    att=list.split(',')
    ind=selAtt+','
    for l in att[:-1]:
        ind=ind+str(l)+','
    ind=ind[:-1]
    indexname = table + '_' + ind
    query = "create index if not exists \"" + indexname + "\" on \"" + table + "\"(" + ind + ");"
    execute_query(conn, query)

def getMVnames(conn):
    return execute_query(conn, "select matviewname from pg_catalog.pg_matviews;")

def dropAllMVs(conn):
    mvnames=getMVnames(conn)
    for n in mvnames:
        ns=[str(i) for i in n]
        execute_query(conn, "drop materialized view \""+ns[0]+"\";")

def getDefOfMV(conn, MVname):
    return execute_query(conn, "select definition from pg_catalog.pg_matviews where matviewname='" + MVname +";")

def getJSONPlannerForQuery(conn, query):
    return execute_query(conn, "explain (format json) " + query)

# generate all cuboids (group bys) including selAtt
def getCuboidsOfAtt(attInGB, selAtt):
    pwset = powerset(attInGB)
    pwset2 = []
    for p in pwset:
        l = list(p)
        l.append(selAtt)
        pwset2.append(tuple(l))
    return pwset2

# percentOfLattice is a float in ]0,1]
# returns nb of created views
# todo for avg, should materialize sum and count
def createMV(conn, attInGB, selAtt, meas, function, table, percentOfLattice,generateIndex=False):
    print("Creating views")
    existing=getMVnames(conn)
    pwset2=getCuboidsOfAtt(attInGB, selAtt)
    # remove last
    del pwset2[-1]

    nbOfMV=len(pwset2)*percentOfLattice
    #print(pwset2)
    #print(int(nbOfMV))

    for i in range(int(nbOfMV)):
        nb = random.randint(0, len(pwset2) - 1)
        gb = pwset2[nb]
        pwset2.remove(gb)
        gbs=''
        for s in gb:
            gbs = gbs + s + ','
        gbs=gbs[:-1]
        #query="create materialized view MV" + str(i) + " as select " + gbs + "," + meas + " from " + table + " group by " + gbs +  ";"
        query="create materialized view \"" + gbs + "\" as select " + gbs + ", " + function + "(" + meas + ") as " + meas + ", count(*) as count  from " + table + " group by " + gbs +  ";"
        #print(query)
        if gbs not in existing:
            execute_query(conn, query)
            if generateIndex==True:
                generateHashIndex(conn, gbs, selAtt)
            if generateIndex=='mc':
                generateMulticolIndex(conn, gbs, gbs, selAtt)
    return nbOfMV

#returns the group by set of a query (having a single group by)
def returnGroupby(query):
    tab=query.split("GROUP BY")
    return tab[1].split(';')[0]




def connect_to_db(dbname, user, password, host='localhost', port='5432'):
    """
    Establishes a connection to the PostgreSQL database.

    :param dbname: Name of the database
    :param user: Database user
    :param password: User's password
    :param host: Database host address (default is 'localhost')
    :param port: Connection port number (default is '5432')
    :return: Connection object
    """
    try:
        conn = psycopg2.connect(
            dbname=dbname,
            user=user,
            password=password,
            host=host,
            port=port
        )
        print("Connection to database established successfully.")
        # Activate extension for sampling rows faster
        act_ext_query = "CREATE EXTENSION  IF NOT EXISTS tsm_system_rows;"
        execute_query(conn, act_ext_query)
        return conn
    except Exception as e:
        print(f"Error connecting to database: {e}")
        return None


def execute_query(conn, query):
    """
    Executes a given SQL query using the established connection.

    :param conn: Connection object
    :param query: SQL query to be executed
    :return: Query result
    """
    try:
        cursor = conn.cursor()
        cursor.execute(query)
        conn.commit()
        #print("Query executed successfully.")

        try:
            result = cursor.fetchall()
            return result
        except psycopg2.ProgrammingError:
            # If the query does not return any data (like an INSERT or UPDATE)
            return None
    except Exception as e:
        print(f"Error executing query: {e}")
        return None
    finally:
        cursor.close()

def printResultSet(result):
    if result is not None:
        for row in result:
            print(row)

def close_connection(conn):
    """
    Closes the database connection.

    :param conn: Connection object
    """
    try:
        conn.close()
        print("Database connection closed.")
    except Exception as e:
        print(f"Error closing connection: {e}")



def generateGB(groupAtts, measures, table):
    """
    Generates group bys and sel from the list of all categorical attributes - unused so far

    :param groupAtts: all the categorical attributes
    :return:
    """
    for g in groupAtts:
        gb = groupAtts.remove(g)
        groupbyAtt = gb
        sel = g
        for m in measures:
            meas = "sum(" + m + ")"

            queryVals = ("select distinct " + sel + " from " + table + ";")


def getSample(conn, measBase, table, sel, sampleSize, method="SYSTEM_ROWS", repeatable=False, valsToSelect=None):
    # sampling using postgresql: https://www.postgresql.org/docs/current/sql-select.html#SQL-FROM
    # system (faster) is block based, bernouili (slower) is row based

    if repeatable:
        is_repeatable = 'REPEATABLE(42)'
    else:
        is_repeatable = ''

    if valsToSelect==None:
        querySample = (
            "SELECT " + sel + ", " + measBase + " FROM " + table + " TABLESAMPLE " + method + " (" + str(
        sampleSize) + ")" + is_repeatable + ";")
    else:
        querySample = (
                "SELECT " + sel + ", " + measBase + " FROM " + table + " TABLESAMPLE " + method + " (" + str(
            sampleSize) + ")" + is_repeatable + " WHERE " + sel + " in " + str(valsToSelect) +";")

    #print('getsample query:', querySample)
    resultVals = execute_query(conn, querySample)
    #print(resultVals)
    #print("Sample size in tuples: ", len(resultVals))
    return resultVals


def emptyGB(conn, nb_of_adom_vals, table, sel, meas):
    #queryEmptyGb = ("SELECT " + sel + ","
    #                + " rank () over (  order by " + meas + " desc ) as rank" +
    #                " FROM " + table + " group by " + sel + " limit " + str(sizeOfVals) + ";")
    #queryEmptyGb = ("SELECT " + sel + ","
    #                + " rank () over (  order by " + meas + " desc ) as rank" +
    #                " FROM " + table + " group by " + sel + ";")
    #queryEmptyGb = ("SELECT " + sel + ","
    #               + " rank () over (  order by " + meas + " desc ) as rank" +
    #              " FROM " + table +  " WHERE " + sel + " in (" + hyp + ") group by " + sel +  ";")

    #hyp = ""
    #for i in range(len(vals)):
    #    hyp = hyp + "'" + str(vals[i]) + "'"
    #    if i != len(vals) - 1:
    #        hyp = hyp + ","

    queryEmptyGb = ("SELECT " + sel + ","
                    + " rank () over (  order by " + meas + " desc ) as rank" +
                    " FROM " + table + " group by " + sel + " limit " + str(nb_of_adom_vals) + ";")

    # No need for this query since we do the complete order next
    #resultEmptyGb = execute_query(conn, queryEmptyGb)

    queryEmptyGbAll = ("SELECT " + sel + ","
                       + " rank () over (  order by " + meas + " desc ) as rank" +
                       " FROM " + table + " group by " + sel + ";")
    #print(queryEmptyGbAll)
    resultEmptyGbAll = execute_query(conn, queryEmptyGbAll)

    return resultEmptyGbAll[:nb_of_adom_vals], resultEmptyGbAll