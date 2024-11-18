import psycopg2
from utilities import powerset
import random
import numpy as np
import pandas as pd
import configparser


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


# returns the set of aggregate queries that can be run over the set of MV in mvnames
def getAggQueriesOverMV(mvnames,sel):
    result=[]
    #print(mvnames)
    for n in mvnames:
        #print("n ",n)
        n=n[0]
        #print("n[0] ",n)
        nset=[]
        for s in n.split(","):
            nset.append(s)
        #print(nset)
        pwset = powerset(nset)
        #print("pwset ",pwset)
        for p in pwset:
            result.append(p)
    #remove duplicates
    #print(result)
    result=list(set(result))
    #remove queries without sel attribute
    #print(result)
    result.remove(())
    #print(result)
    res=[]
    for r in result:
        #print(r)
        #print(r[-1])
        if r[-1]==sel:
            res.append(r)
    #print(res)
    return res

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
                #print('creating index on view')
                generateHashIndex(conn, gbs, selAtt)
            if generateIndex=='mc':
                #print('creating multicolumn index on view')
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


def generateArtificialDataset_V1(conn,num_rows = 50000, nbAtt=10):
    # Generating values for the first attribute (categorical, repeated, exponential distribution)
    exp_values = np.random.exponential(scale=1.0, size=num_rows)
    first_attribute = pd.qcut(exp_values, q=10, labels=[f'Category_{i+1}' for i in range(10)])

    # Generating values for the second attribute (numerical, Gaussian distribution)
    gaussian_values = np.random.normal(loc=50, scale=10, size=num_rows)

    # Generating values for the other attributes (categorical, repeated, uniformly random)
    categorical_values = [f'Category_{i+1}' for i in range(15)]
    other_attributes = {
        f'A_{i+3}': np.random.choice(categorical_values, size=num_rows)
        for i in range(nbAtt)
    }

    # Creating the DataFrame
    data = {
        'Attribute_1': first_attribute,
        'Attribute_2': gaussian_values,
    }
    data.update(other_attributes)

    # Creating the DataFrame and saving to CSV
    df = pd.DataFrame(data)
    path='/Users/marcel/PycharmProjects/hoeffding4ranks/data/'
    df.to_csv(path+'relational_table.csv', index=False)

    print("Relational table with ",num_rows," rows and ", nbAtt+2," attributes has been generated and saved as 'relational_table.csv'")

    tablename="artificial_"+str(nbAtt)
    query = "drop table if exists \"" + tablename + "\";"
    execute_query(conn, query)
    att='Attribute_1 varchar, Attribute_2 float, '
    for i in range(1,nbAtt+1):
        att=att+"A_"+str(i+2)+" varchar, "
    att=att[:-2]
    query="create table "+tablename+"("+att+");"
    execute_query(conn, query)
    query = "copy " + tablename + " from \'" + path+'relational_table.csv' + "\' (header, format csv);"
    print(query)
    execute_query(conn, query)



def generateArtificialDataset(conn,num_rows = 50000, nbAtt=10):

    # Set the number of unique categories for the categorical attributes
    num_categories_first_attr = 10
    num_categories_other_attrs = 5

    # 1st Attribute: Categorical values with occurrences drawn from an exponential distribution
    categories_first_attr = [f'Category_{i+1}' for i in range(num_categories_first_attr)]
    occurrences = np.random.exponential(scale=num_rows / num_categories_first_attr, size=num_categories_first_attr).astype(int)
    occurrences = occurrences * num_rows // occurrences.sum()  # Scale occurrences to match the total number of rows

    first_attribute = []
    for i, count in enumerate(occurrences):
        first_attribute.extend([categories_first_attr[i]] * count)
    #this is because GPT does not know how to count
    while len( first_attribute) != num_rows:
        first_attribute.append('Category_1')
    # Shuffle to randomize order
    random.shuffle(first_attribute)

    # 2nd Attribute: Numerical values drawn from a Gaussian (normal) distribution
    mean = 50
    std_dev = 10
    second_attribute = np.random.normal(loc=mean, scale=std_dev, size=num_rows)

    # Other Attributes: Categorical values drawn uniformly at random
    other_attributes = []
    for i in range(nbAtt):  # nbAtt other categorical attributes
        other_attribute = [f'Cat_{i+1}_{j+1}' for j in range(num_categories_other_attrs)]
        other_attributes.append(np.random.choice(other_attribute, size=num_rows))

    # Create a DataFrame to represent the relational table
    data = {
        'Attribute_1': first_attribute,
        'Attribute_2': second_attribute
    }

    for i in range(nbAtt):
        data[f'A_{i+3}'] = other_attributes[i]

    #print(data)
    df = pd.DataFrame(data)

    # Display the DataFrame structure
    df.info()

    # Optionally, save the DataFrame to a CSV file
    #df.to_csv('generated_table_instance.csv', index=False)


    df = pd.DataFrame(data)
    path='/Users/marcel/PycharmProjects/hoeffding4ranks/data/'
    df.to_csv(path+'relational_table.csv', index=False)

    print("Relational table with ",num_rows," rows and ", nbAtt+2," attributes has been generated and saved as 'relational_table.csv'")

    tablename="artificial_"+str(nbAtt)
    query = "drop table if exists \"" + tablename + "\";"
    execute_query(conn, query)
    att='Attribute_1 varchar, Attribute_2 float, '
    for i in range(1,nbAtt+1):
        att=att+"A_"+str(i+2)+" varchar, "
    att=att[:-2]
    query="create table "+tablename+"("+att+");"
    execute_query(conn, query)
    query = "copy " + tablename + " from \'" + path+'relational_table.csv' + "\' (header, format csv);"
    #print(query)
    execute_query(conn, query)


if __name__ == "__main__":

    config = configparser.ConfigParser()

    # The DB wee want
    #config.read('configs/flights1923.ini')
    #config.read('configs/flights.ini')
    config.read('configs/artificial.ini')
    #config.read('configs/ssb.ini')
    # The system this is running on
    USER = "PM"

    # Database connection parameters
    dbname = config[USER]['dbname']
    user = config[USER]['user']
    password = config[USER]['password']
    host = config[USER]['host']
    port = int(config[USER]['port'])


    # Connect to the database
    conn = connect_to_db(dbname, user, password, host, port)

    dropAllMVs(conn)
    generateArtificialDataset(conn)